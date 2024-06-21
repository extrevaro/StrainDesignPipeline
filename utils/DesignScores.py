import time

import itertools
import pandas as pd
import numpy as np
import cobra as cb
#import multiprocessing as mp
import joblib

from cobra.flux_analysis import production_envelope, flux_variability_analysis
from cameo.flux_analysis.simulation import pfba

def get_metabolic_intersections(rxn_list):
    met_dict = dict()
    
    for r1, r2 in itertools.combinations(rxn_list, 2):
        intersecting_mets = set([p.id for p in r1.products]).intersection(set([s.id for s in r2.reactants]))
        met_dict['|'.join([r1.id, r2.id])] = intersecting_mets
        
    return met_dict

def get_total_carbon_uptake(model):

    carbon_uptakes = [ model.metabolites.get_by_id(m).elements['C']*u
                       if 'C' in model.metabolites.get_by_id(m).elements.keys()
                       else 0
                       for m,u in zip(model.summary().uptake_flux['metabolite'], model.summary().uptake_flux['flux'])
                     ]
    
    return sum(carbon_uptakes)

def pathway_yield(model_list, pathway_rxn_list):
    #define a model that will consist on the fusion of the active metabolic network
    #of each model contributing in the distribution of the target pathway
    mixed_model = cb.Model('mixed_model')
    #get pathway reaction conections through intermediate metabolites
    pathway_intersections = get_metabolic_intersections(pathway_rxn_list)
    #get those metabolites connecting the pathway
    conected_metabolites = set([m for m_s in pathway_intersections.values() for m in m_s])
    
    for model in model_list:
        #find target reactions and metabolites within the present model
        targets_in_model = [model.reactions.get_by_id(r.id)
                            for r in pathway_rxn_list
                            if r.id in [rxn.id for rxn in model.reactions]]
        
        mets_in_model = set([m.id
                             for r in targets_in_model
                             for m in r.metabolites
                             if m.id in conected_metabolites])
        #get the metabolites not needing imports for the model to produce them 
        #(are conected with the metabolic network)
        intersections_in_model = get_metabolic_intersections(targets_in_model)
        conected_metabolites_in_model = set([m for m_s in intersections_in_model.values() for m in m_s])
        
        #include exchange reaction to disconected metabolites to nable synthesis
        #by the present model
        for met in mets_in_model-conected_metabolites_in_model:
            #only follow the flux of carbon-containing metabolites
            if 'C' in model.metabolites.get_by_id(met).elements.keys():
                exchange_out = cb.Reaction('_'.join(['EX', met, 'out_artificial']))
                exchange_out.add_metabolites({model.metabolites.get_by_id(met): -1.0})
                exchange_in = cb.Reaction('_'.join(['EX', met, 'in_artificial']))
                exchange_in.add_metabolites({model.metabolites.get_by_id(met): 1.0})
                model.add_reactions([exchange_in, exchange_out])           

        #force the production of all intermediates in the model
        for r in targets_in_model:
            model.reactions.get_by_id(r.id).bounds = (0.05, 1000)
            
        try:
            #get all active reactions for this model configuration through a fva
            #model.objective = targets_in_model[-1]
            fva_result = flux_variability_analysis(model, model.reactions, fraction_of_optimum=0.5, processes=1)
            active_reactions = fva_result.loc[(fva_result['minimum']!=0)|(fva_result['maximum']!=0)].index.tolist()
            #fba_result = model.optimize().to_frame()
            #active_reactions = fba_result.loc[fba_result['fluxes']!=0].index.tolist()
            #construct the reaction network without the added exchanges of disconected metabolites
            #as they need to be in the other model for the pathway to work, if not the pathway
            #distribution will result in an unfeasible probustion of the target
            new_reactions = [model.reactions.get_by_id(r)
                             for r in active_reactions
                             if 'artificial' not in r
                             and r not in [rxn.id for rxn in mixed_model.reactions]]
            
            mixed_model.add_reactions(new_reactions)
            
            for r in mixed_model.reactions:
                if r.id in [target.id for target in targets_in_model]:
                    r.bounds = (fva_result.loc[r.id].minimum, fva_result.loc[r.id].maximum)
            
        except Exception as e:
            print(e)
            result = 0
            
    #get the maximum yield for this pathway distribution using the mixed model
    mixed_model.objective = pathway_rxn_list[-1].id
    target_metabolite = [m for m in pathway_rxn_list[-1].reactants][0]
    pathway_yield = mixed_model.slim_optimize()
    carbon_product_flux = target_metabolite.elements['C']*pathway_yield
    
    #to measure efficiency of the pathway distribution we have to scale the
    #maximum yield with the initial carbon uptake
    carbon_uptake = get_total_carbon_uptake(mixed_model)
    carbon_yield = carbon_product_flux/carbon_uptake
            
    return carbon_yield

def coupling_strength(model, target_reaction, target_biomass, carbon_source):

    prod_env = production_envelope(model, 
                                   [target_biomass], 
                                   objective=target_reaction, 
                                   carbon_sources=carbon_source)
    
    '''
    AS the score is compute differently depending if reaction is coupled, first compute
    the decision parameter to check if solution couple or not the production with the
    growth. Coupling means that in the envelope, the minimum of the points where the 
    biomass is equal to its maximum should be above 0 (there is production at maximum
    growth).This value is computed in the following variable:
    '''
    max_growth = max(prod_env[target_biomass].tolist())
    min_p_max_growth = min(prod_env.loc[prod_env[target_biomass]==max_growth,
                                                 'flux_minimum' ].tolist())
    
    is_coupled = min_p_max_growth > 0
    
    
    if is_coupled:
        '''
        Before the computation compute the rest of needed parameters:
        min_p_no_growth : minimum production of target when there is no growth
        max_growth_no_p : maximum growth when there is no production
        '''
        min_p_no_growth = min(prod_env.loc[prod_env[target_biomass]==0,
                                                    'flux_minimum' ].tolist())
        
        min_p = min(prod_env['flux_minimum'].tolist())
        if min_p > 0:
            max_growth_no_p = 0
        else:
            max_growth_no_p = max(prod_env.loc[prod_env['flux_minimum']==0,
                                               target_biomass].tolist())
            
        coupling_strength = min([1,(min_p_no_growth/min_p_max_growth)]) + ((max_growth-max_growth_no_p)/max_growth)
        
    if not is_coupled:
        '''
        Before the computation compute the rest of needed parameters:
        max_growth_delta_p : maximum growth when the production is forced to a
                             minimum value above 0, called delta
        
        Here delta will be the minimum nonzero value in the flux_maximum column
        '''
        #check if max value of 'flux_maximum' is higher than 0 to avoid error
        if prod_env['flux_maximum'].max()>0:
            delta = min(prod_env.loc[prod_env['flux_maximum']>0, 'flux_maximum'].tolist())
            max_growth_delta_p = max(prod_env.loc[prod_env['flux_maximum']==delta,
                                                          target_biomass].tolist())
                                                                  
            coupling_strength = min([0, (max_growth-max_growth_delta_p)/delta])
            
        else:
            coupling_strength = 0
    
    return coupling_strength, is_coupled
    
class StrainDesignScorer:
    def __init__(self, model_setup, sd, score, parallelize):
        self.model_setup = model_setup
        self.sd = sd
        self.score = score
        self.parallelize = parallelize
        
        
    def parallel_pfba_yield_analysis(self, iterable):
        model=self.model_setup['model']
        target_carbons = [m for m in model.reactions.get_by_id(self.model_setup['target_reaction']).metabolites][0].elements['C']
        uptake_carbons = [m for m in model.reactions.get_by_id(self.model_setup['carbon_source']).metabolites][0].elements['C']
        duration = 0

        for rxn, bounds in self.model_setup['model_bounds'].items():
            model.reactions.get_by_id(rxn).bounds = bounds

        for d, r in zip(iterable==1, np.array(self.model_setup['reaction_list'])):
            if d:
                model.reactions.get_by_id(r).knock_out()
        try:
            pfba_result = pfba(model)

            if pfba_result[self.model_setup['target_biomass']] < 0.05*self.model_setup['max_target_biomass_flux']: #only compute yield if growt is above 5 percent of WT
                product_yield = 0
            
            else:
                target_flux = pfba_result[self.model_setup['target_reaction']]            
                uptake_flux = abs(pfba_result[self.model_setup['carbon_source']])                    
                product_yield = (target_flux*target_carbons)/(uptake_flux*uptake_carbons)

            if product_yield > 0 :
                end = time.time()
                duration = end-self.start           

        except Exception as e:
            product_yield = 0
            
        yield_result =list(iterable)+[product_yield]
        #print(yield_result)

        return yield_result, duration
        
    
    def parallel_pfba_bcpy_analysis(self, iterable):
        model=self.model_setup['model']
        duration = 0

        for rxn, bounds in self.model_setup['model_bounds'].items():
            model.reactions.get_by_id(rxn).bounds = bounds

        for d, r in zip(iterable==1, np.array(self.model_setup['reaction_list'])):
            if d:
                model.reactions.get_by_id(r).knock_out()
        try:
            pfba_result = pfba(model)
            bcpy = pfba_result[self.model_setup['target_reaction']]*pfba_result[self.model_setup['target_biomass']]

            if bcpy > 0 :
                end = time.time()
                duration = end-self.start

        except Exception as e:
            print(e)
            bcpy = 0
            
        bcpy_result =list(iterable)+[bcpy]

        return bcpy_result, duration
    
    
    def parallel_coupling_strength(self, iterable):
        model = self.model_setup['model']
        duration = 0

        for rxn, bounds in self.model_setup['model_bounds'].items():
            model.reactions.get_by_id(rxn).bounds = bounds

        for d, r in zip(iterable==1, np.array(self.model_setup['reaction_list'])):
            if d:
                model.reactions.get_by_id(r).knock_out()
        try:
            c_s = coupling_strength(model,
                                    self.model_setup['target_reaction'],
                                    self.model_setup['target_biomass'],
                                    self.model_setup['carbon_source'])[0]

            if c_s > 0 :
                end = time.time()
                duration = end-self.start

        except Exception as e:
            c_s = 0
            
        cs_result = list(iterable)+[c_s]

        return cs_result, duration
    
    
    def parallel_gcfront_score(self, iterable):
        model = self.model_setup['model']
        duration = 0
        max_c_s = 2

        for rxn, bounds in self.model_setup['model_bounds'].items():
            model.reactions.get_by_id(rxn).bounds = bounds

        for d, r in zip(iterable==1, np.array(self.model_setup['reaction_list'])):
            if d:
                model.reactions.get_by_id(r).knock_out()
        try:
            pfba_result = pfba(model)

            target_flux = pfba_result[self.model_setup['target_reaction']]
            target_biomass_flux = pfba_result[self.model_setup['target_biomass']]
            growth_above_threshold = target_biomass_flux > (0.1*self.model_setup['max_target_biomass_flux'])
            c_s, is_coupled = coupling_strength(model,
                                                self.model_setup['target_reaction'],
                                                self.model_setup['target_biomass'],
                                                self.model_setup['carbon_source'])

            if is_coupled :
                end = time.time()
                duration = end-self.start
                
            gcfront_score = (target_flux/self.model_setup['max_target_flux'])*int(growth_above_threshold) + \
            		     (target_biomass_flux/self.model_setup['max_target_biomass_flux'])*int(is_coupled) + \
            		     (c_s/max_c_s)*int(growth_above_threshold)

        except Exception as e:
            target_flux = 0
            target_biomass_flux = 0
            c_s = 0
            is_coupled = 0
            gcfront_score = 0
            
        
        gcfront_result = list(iterable)+[gcfront_score]

        return gcfront_result, duration
    
    
    def parallel_gcfront_components(self, iterable):
        model = self.model_setup['model']
        duration = 0

        for rxn, bounds in self.model_setup['model_bounds'].items():
            model.reactions.get_by_id(rxn).bounds = bounds

        for d, r in zip(iterable==1, np.array(self.model_setup['reaction_list'])):
            if d:
                model.reactions.get_by_id(r).knock_out()
        pfba_result = pfba(model)

        target_flux = pfba_result[self.model_setup['target_reaction']]
        target_biomass_flux = pfba_result[self.model_setup['target_biomass']]        
        
        try:
            c_s, is_coupled = coupling_strength(model,
                                                self.model_setup['target_reaction'],
                                                self.model_setup['target_biomass'],
                                                self.model_setup['carbon_source'])

            if is_coupled :
                end = time.time()
                duration = end-self.start

        except Exception as e:
            c_s = 0
            is_coupled = 0

        gcfront_result = (target_biomass_flux, target_flux, c_s)

        return gcfront_result, duration
        

    def parallel_pathway_distribution_yield(self, iterable):
    	assert type(self.model_setup['model']) == list and len(self.model_setup['model'])>1,"For distributing a pathway you need to pass a list of models."
    	
    	assert all([type(r)==cb.core.reaction.Reaction for r in self.model_setup['reaction_list']]),"Target list need to be composed of cobra.Reaction elements!"
    	
    	assert len(iterable)==len(self.model_setup['model'])*len(self.model_setup['reaction_list']),"Iterable need to be composed of (n_models x n_targets) elements"
    	
    	model_index = 0
    	iterable_pos = 0
    	rxns_to_add = []
    	model_list = []

    	for d, reaction in zip(iterable, np.array(self.model_setup['reaction_list']*len(self.model_setup['model']))):
            
            if d==1:
                rxns_to_add.append(reaction)
            
            iterable_pos += 1

            if iterable_pos % len(self.model_setup['reaction_list']) == 0:
                model = cb.Model('MODEL_%u' % model_index)
                model_reactions = [r for r in self.model_setup['model'][model_index].reactions]
                rxns_to_add += model_reactions
                model.add_reactions(rxns_to_add)
                model_list.append(model)
                del model               
                model_index += 1
                rxns_to_add = []         
        
        #Now execute a function that taking model_list computes the pathway yield
    	try:
            distribution_yield = pathway_yield(model_list, self.model_setup['reaction_list'])            

    	except Exception as e:
            print(e)
            distribution_yield = 0

    	return distribution_yield
    
    
    def run(self):
        self.start = time.time()

        if self.parallelize:
            #in order to work in notebook
            '''
            try:
                mp.set_start_method("spawn", force=True)
        
            except RuntimeError:
                pass
            '''
            #parallelize work
            print('Testing all samples in the GEM...')
            print('Parallelizing task..')
            processes = 6
            #pool = mp.Pool(processes)
            
            if self.score == 'bcpy':
                score_model = self.parallel_pfba_bcpy_analysis

            elif self.score == 'coupling_strength':
                score_model = self.parallel_coupling_strength

            elif self.score == 'gcfront':
                score_model = self.parallel_gcfront_score
                
            elif self.score == 'components':            
                score_model = self.parallel_gcfront_components

            elif self.score == 'pathway_distribution':            
                score_model = self.parallel_pathway_distribution_yield

            elif self.score == 'yield':                
                score_model = self.parallel_pfba_yield_analysis
                
            else:
                raise ValueError('Invalid score parameter')
                print('Options are : yield, bcpy, coupling_strength, gcfront or components')

            try:
                result_list = joblib.Parallel(backend="loky", n_jobs=processes, batch_size=round(len(self.sd)/processes), verbose=10)(joblib.delayed(score_model)(r) for r in self.sd)
                #result_list = pool.map(score_model, self.sd)
                #pool.close()
                #pool.join()

            except Exception as e:
                print(e)
                #pool.close()
                #pool.join()
                print('Parallelizing error!!')
                
            time_list=np.array([r[1] for r in result_list])
            gc_time_list = time_list[time_list>0]

            print('End of parallel task!')
            print('A total of %u Growth-Copupled designs were found!' % len(gc_time_list))
            
            if len(gc_time_list)>0:
                print('First design found in %f seconds' % min(gc_time_list))
            
            return result_list
            
            
        else:
        
            if self.score == 'bcpy':
                return self.parallel_pfba_bcpy_analysis(self.sd)

            elif self.score == 'coupling_strength':
                return self.parallel_coupling_strength(self.sd)

            elif self.score == 'gcfront':
                return self.parallel_gcfront_score(self.sd)
            
            elif self.score == 'components':
                return self.parallel_gcfront_components(self.sd)
                
            elif self.score == 'pathway_distribution':            
                return self.parallel_pathway_distribution_yield(self.sd)
                
            elif self.score == 'yield':                
                return self.parallel_pfba_yield_analysis(self.sd)
                
            else:
                raise ValueError('Invalid score parameter')
                print('Options are : yield, bcpy, coupling_strength, gcfront or components')
                
                
class objective_function:
    def __init__(self, model_setup, target_score, max_dels):
        self.model_setup=model_setup
        self.target_score=target_score
        self.max_dels=max_dels
    
    def compute(self, Designs:np.array):
        result = []
        # Your binary function here
        for design in Designs:
            if design.sum() <= self.max_dels:
                scorer = StrainDesignScorer(self.model_setup, design,
                                            score=self.target_score, parallelize=False)

                result.append(-scorer.run()[0][-1])

            else:
                result.append(0.0)

        return result

