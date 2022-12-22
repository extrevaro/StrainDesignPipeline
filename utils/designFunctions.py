from cameo.flux_analysis.simulation import lmoma, pfba
from cameo.strain_design.deterministic.linear_programming import OptKnock
from cameo.flux_analysis.analysis import phenotypic_phase_plane
#FSEOF ANALYSIS
from cameo.visualization.plotting.with_plotly import PlotlyPlotter
from cameo.strain_design.deterministic.flux_variability_based import FSEOF
import cobra
from cobra.sampling import OptGPSampler, ACHRSampler
from cobra import Reaction
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis.variability import find_essential_genes
import pandas as pd
import numpy as np
import straindesign as sd
import plotly.express as px
import plotly.graph_objects as go
import ast
import re
import os
import sys
import time
import random

def get_rxn_with_fva_flux(model, fraction_of_objective=0.2, loopless=False):
    model = model.copy()
    #get fva results for all reactions in the model
    fva_result = cobra.flux_analysis.variability.flux_variability_analysis(model, fraction_of_optimum=fraction_of_objective,  loopless=loopless)
    #Filter out those reaction tath have min and max fluxes equals to 0:
    rxn_with_flux = fva_result.loc[(fva_result['maximum'] != 0) | (fva_result['minimum'] != 0) ].index.tolist()
    return rxn_with_flux

def custom_minimal_media(model, carbon_sources):
    from cobra.medium import minimal_medium
    max_growth = model.slim_optimize()
    minimal_media = minimal_medium(model, max_growth)
    custom_media = { r: (minimal_media[r], model.reactions.get_by_id(r).upper_bound) for r in minimal_media.index }
    for cs in carbon_sources.keys():
        custom_media[cs] = (carbon_sources[cs], model.reactions.get_by_id(cs).upper_bound)
    
    model.medium = {r: abs(custom_media[r][0]) for r in custom_media.keys()}
    
    return model, custom_media

def display_KO_candidates_results(model, ko_targets, target_biomass, target_reaction, engine, media=None):
    for s in ko_targets:
        print('Results for deletion of %s' % ', '.join([r for r in s]))
        with model:
            for r in s:
                model.reactions.get_by_id(r).knock_out()
            if engine == 'cameo':
                from cameo import phenotypic_phase_plane
                from cameo.visualization.plotting.with_plotly import PlotlyPlotter
                plotter = PlotlyPlotter()
                pfba_solution = pfba(model)
                #lmoma_solution = lmoma(model, reference=pfba_solution.fluxes)
                mutant_biomass = pfba_solution.fluxes[target_biomass]
                mutant_target = pfba_solution.fluxes[target_reaction]
                print('Flux of biomass solution is %s and the flux of target reaction is %s' % (mutant_biomass, mutant_target))
                try:
                    result = phenotypic_phase_plane(model,
                                                    variables=[model.reactions.get_by_id(target_biomass)],
                                                    objective=model.reactions.get_by_id(target_reaction))
                    result.plot(plotter)
                except:
                    print('Infeasible solution, trying_next')
                    continue

            elif engine == 'mewpy':
                from mewpy.visualization.envelope import plot_flux_envelope
                from mewpy.simulation import get_simulator
                if media:
                    simul = get_simulator(model, envcond=media)
                    constraints  = { r :  (0,0) for r in s}
                    lmoma_solution = simul.simulate(method=SimulationMethod.lMOMA, constraints=constraints)
                    mutant_biomass = lmoma_solution.fluxes[target_biomass]
                    mutant_target = lmoma_solution.fluxes[target_reaction]
                    print('Flux of biomass solution is %s and the flux of target reaction is %s' % (mutant_biomass, mutant_target))
                    try:
                        display(plot_flux_envelope(simul, target_biomass, target_reaction, constraints = constraints))
                    except:
                        print('Infeasible solution, trying_next')
                        continue
                else:
                    print('Production envelope could not be plotted in mewpy without the media specification!')

            else:
                print('You have to select between cameo and mewpy for variable engine if you want to see the production envelope!')

            for r in s:
                display(model.reactions.get_by_id(r))

def get_expression_modulation_targets(result_df, threshold, n_of_results=10):
    upper_threshold = threshold
    lower_threshold = 1/threshold
    downregulated_genes = [index for index, row in result_df.iterrows() if (row[n_of_results]/row[1])<lower_threshold]
    upregulated_genes = [index for index, row in result_df.iterrows() if (row[n_of_results]/row[1])>upper_threshold]
    return downregulated_genes, upregulated_genes

def display_unique_strategies(results_list):
    unique_strategies = []
    for result in results_list:
        for strategy in result.data_frame.reactions:
            if strategy not in unique_strategies:
                unique_strategies.append(strategy)
    return unique_strategies


#this function will create the protected reactions.
#This set consists on those reactions whose deletion will seriously hamper the flux
#of both target biomass and metabolite exchange
def select_protected_reactions(model, target_biomass, target_reaction, reactions_to_test=None, flux_fraction=0.1):
    wt_sol = model.optimize()
    growth_threshold = wt_sol.objective_value*flux_fraction
    minimum_production = wt_sol.fluxes[target_reaction]*flux_fraction
    protected_reactions = []

    if reactions_to_test:
        deletion_list = [r for r in model.reactions if r.id in reactions_to_test]
    else:
        deletion_list = model.reactions

    print('Checking %s reactions...' % len(deletion_list))
    for r in deletion_list:
        test_model = model.copy()
        model.objective = target_biomass
        test_model.reactions.get_by_id(r.id).knock_out()
        try:
            pfba_result = pfba(test_model)
            no_flux_biomass = pfba_result[target_biomass] < growth_threshold
            no_flux_metabolite = pfba_result[target_reaction] < minimum_production

        except Exception as e:
            print(e)
            print(r.id)
            protected_reactions.append(r.id)
            continue

        if no_flux_biomass or no_flux_metabolite:
            print(r.id)
            protected_reactions.append(r.id)

    return set(protected_reactions)

def delete_trasport_reactions(model, rxn_list):
    filtered_list=[]
    for r in rxn_list:
        sustrate_compartiments = {met.id[:-2] : met.id[-2:] for met in model.reactions.get_by_id(r).reactants}
        product_compartiments = {met.id[:-2] : met.id[-2:] for met in model.reactions.get_by_id(r).products}
        shared_mets = set(sustrate_compartiments.keys()) & set(product_compartiments.keys())
        if all([sustrate_compartiments[met] == product_compartiments[met] for met in shared_mets]):
            filtered_list.append(r)
    print(len(filtered_list))
    return filtered_list

def find_coupled_reactions(model, reaction_to_test):
    """Find reaction sets that are structurally forced to carry equal flux"""
    model = model.copy()
    stoichiometries = {}
    for reaction in reaction_to_test:
        for met, coef in model.reactions.get_by_id(reaction).metabolites.items():
            stoichiometries.setdefault(met.id, {})[reaction] = coef

    # Find reaction pairs that are constrained to carry equal flux
    couples = []
    for met_id, stoichiometry in stoichiometries.items():
        if len(stoichiometry) == 2 and set(stoichiometry.values()) == {1, -1}:
            couples.append(set(stoichiometry.keys()))

    # Aggregate the pairs into groups
    coupled_groups = []
    for couple in couples:
        for group in coupled_groups:
            if len(couple & group) != 0:
                group.update(couple)
                break
        else:
            coupled_groups.append(couple)

    return coupled_groups

#function for generate target list of knocks outs.
#Uses a model and a reaction set of protected reactions computed by the function
#select_protected_reactions. For coherent results when used with no blocked_reactions
#(default), it should be called with the same epsilon and tolerance used to generate
#the consistent model through fastcc
def get_KO_candidate_list(model, protected_reactions, blocked_reactions=None, carbon_limit=16, epsilon=1e-4, tolerance=1e-7):
    #exclude protected_reactions:
    candidate_set = set([r.id for r in model.reactions]) - protected_reactions
    #exclude blocked reactions
    if not blocked_reactions:
        blocked_reactions = Fastcore.fast_find_blocked(model, epsilon=epsilon, tolerance=tolerance)
    candidate_set = candidate_set - blocked_reactions
    #exclude non GPR, exchange and spontaneous reactions
 
    candidate_list = [ r for r in candidate_set if
                      model.reactions.get_by_id(r).gene_reaction_rule != '' and #Non GPR
                      len(model.reactions.get_by_id(r).metabolites) > 1 and     #Exchange
                      'spontaneous' not in model.reactions.get_by_id(r).name ]  #Spontaneous 
    
    #exclude transport reactions
    candidate_list = delete_trasport_reactions(model, candidate_list)
    #exclude reactions using high carbon molecules (>16 C atoms)
    
    cofactor_list = ['atp', 'adp', 'coa', 'nadph', 'nadp', 'nadh', 'nad', 'fad', 'fadh2']
    high_carbon_reactions=[]
    for r in candidate_list:
        for subs in model.reactions.get_by_id(r).reactants:
            formula = subs.formula
            id = subs.id
            if not any([cofactor in id for cofactor in cofactor_list]):
                try:
                    n_of_carbons = int(re.findall(r'C([0-9]{1,2})', formula)[0])
                except:
                    continue
                if n_of_carbons > carbon_limit:
                    high_carbon_reactions.append(r)
                    break
    

    candidate_list = list(set(candidate_list)-set(high_carbon_reactions))
    

    #delete reactions that are coupled with others in the list
    
    rxns_groups = find_coupled_reactions(model, candidate_list)
    candidate_list = [ r for r in candidate_list if
                       r not in [not_leader for rxn_group in rxns_groups for not_leader in list(rxn_group)[1:]] ]

    return candidate_list

def execute_OptKnock(configured_model, m_ko, target_reaction, target_biomass,to_exclude):
    optknock = OptKnock(configured_model, exclude_reactions=to_exclude,
                        remove_blocked=False,
                        exclude_non_gene_reactions=False,
                        fraction_of_optimum=0.1,
                        use_nullspace_simplification=False)
    
    optknock_result = optknock.run(max_knockouts=m_ko, target=target_reaction, biomass=target_biomass, max_results=10)
    return optknock_result


def generate_strategies(set_up_params, configured_model, target_biomass, target_reaction, carbon_source, blocked_reactions):
    #compute reference values for biomas and target fluxes:
    experiment = '_'.join(set_up_params['project'].split('_')[0:-2])   
    #try to read file containing essential reactions, if not possible compute them
    essential_gene_filename = experiment+'_essential_reactions.csv'
    project_path = '/'.join([set_up_params['dir_path'], set_up_params['project']])
    framework_path = '/'.join([project_path, set_up_params['framework']])
    essential_gene_filepath = '/'.join([project_path, essential_gene_filename])
    try:
        essential_df = pd.read_csv(essential_gene_filepath, header=None)
    except IOError:
        print("File not accessible, computing process esential reactions...")
        essential_rxns = select_protected_reactions(configured_model, target_biomass, target_reaction)
        pd.DataFrame.from_dict({'Reaction_id': list(essential_rxns)}).to_csv(essential_gene_filepath, index=False)
        #make directory for the project inside the current directory:
        if not os.path.isdir(project_path):
            os.mkdir(project_path)
                            
        essential_df = pd.read_csv(essential_gene_filepath, header=None)

    essential_rxns = set([r[0] for r in essential_df.values.tolist()])
    
    print('Essential reactions for bioprocess are: \n %s' % (essential_rxns))
    #make directory for the framework inside the current directory:
    if not os.path.isdir(framework_path):
        os.mkdir(framework_path)
    #get the candidate list:                          
    for c_l in set_up_params['max_cl_range']:
        for m_ko in set_up_params['max_knock_out_range']:
            print(c_l)
            candidate_rxns = get_KO_candidate_list(configured_model, essential_rxns, blocked_reactions=blocked_reactions, carbon_limit=c_l)
            
            print('All candidate reactions are :')
            print(candidate_rxns)

            for r in range(set_up_params['replicates']):
                filename= '_'.join([experiment, str(m_ko), 'ko', str(c_l),str(r),'rep.csv'])
                filepath= '/'.join([framework_path, filename])
                candidate_subset = set([ candidate_rxns[random.randint(0,len(candidate_rxns)-1)] for i in range(0, 130) ])
                print('Searching for Knock Out strategies using a list of %s candidates' % len(candidate_subset))
                if set_up_params['framework'] == 'cameo':
                    to_exclude = set([r.id for r in configured_model.reactions])- candidate_subset
                    OptKnock_result = execute_OptKnock(configured_model, m_ko, target_reaction, target_biomass,to_exclude)
                    OptKnock_result.data_frame.to_csv(filepath)
                    

def optknock_result_analysis(configured_model, target_biomass, carbon_source, target_reaction,
                             directory_name, max_cl_range=range(12,6,-1), max_knock_out_range=range(10,16,1),
                             top_results = 20, performance_parameter='carbon_yield_maximum'):
    #get the main features of results and compute carbon yield for every strategy found in OptKnock:
    cl_list = []
    m_ko_list = []
    cy_list = []
    mb_list =[]
    strategy_list = []
    organism = directory_name.split('_')[0]
    for c_l in max_cl_range:
        for m_ko in max_knock_out_range:
            filename= organism+'_'+str(m_ko)+'_ko_'+str(c_l)+'_'+str(r)+'_rep.csv'
            filepath= '/'.join([directory_name, filename])
            print('Reading %s' % filepath)
            df = pd.read_csv(filepath)

            if len(df) > top_results:
                df = df.sort_values([target_reaction, 'biomass'], ascending = [False, True])[:top_results]

            for s in df.reactions.values.tolist():
                test_model = configured_model.copy()
                strategy = list(ast.literal_eval(s))
                test_model.remove_reactions(strategy)
                prod_env = production_envelope(test_model, [target_biomass], objective=target_reaction, carbon_sources=carbon_source)
                #get the carbon yield and biomass and save to a variable to later plot them
                max_biomass = max(prod_env[target_biomass].values.tolist())
                carbon_yield = prod_env.loc[prod_env[target_biomass] == max_biomass, performance_parameter].values[0]
                print(carbon_yield)
                strategy_list.append(strategy)
                cl_list.append(c_l)
                m_ko_list.append(m_ko)
                cy_list.append(carbon_yield)
                mb_list.append(max_biomass)
                
    data = {'Strategy' : strategy_list,
            'Max_KO' : m_ko_list,
            'Carbon_limit' : cl_list, 
            'Max_biomass' : mb_list, 
            'Max_carbon_yield' : cy_list}

    results_df = pd.DataFrame.from_dict(data)
    return results_df



def wait_for_file(fn: str, max_wait_sec: int = 3600 * 24 * 365,
                    check_interval: float = 0.02) -> bool:
    """
    Waits for file maximum of max_wait_sec. Returns True if file was detected within specified max_wait_sec
    Args:
      fn: filename on task machine
      max_wait_sec: how long to wait in seconds
      check_interval: how often to check in seconds
    Returns:
      False if waiting was was cut short by max_wait_sec limit, True otherwise
    """
    print("Waiting for file", fn)
    start_time = time.time()
    while True:
        if time.time() - start_time > max_wait_sec:
            return False
        if not os.path.exists(fn):
            time.sleep(check_interval)
            continue
        else:
            break
    return True 



def flux_sampling_screening(comb_deletions, df, model, target_biomass, target_reaction, flux_fraction=0.2):
    wt_growth = float(model.optimize().objective_value)
    growth_threshold = wt_growth*flux_fraction
    print('Screening a total of %s strategies using flux sampling, this is going to take a while...' % str(len(comb_deletions)))
    for del_comb in comb_deletions:
        test_model = model.copy()
        mut_name = '-'.join(del_comb)
        for deletion in del_comb:
            test_model.reactions.get_by_id(deletion).knock_out()
            
        #To ensure flux sampling results have a minimum biomass of 20% of wt,
        #change the lower bound to this value
        test_model.reactions.get_by_id(target_biomass).bounds = (growth_threshold, 1000)
        try:
            optgp = OptGPSampler(test_model, processes=8, thinning=1000)
            data = optgp.sample(5600)
            data_dict = {target_biomass : data[target_biomass], 
                         target_reaction : data[target_reaction],
                         'strain' : [mut_name]*5600}        
            comb_df = pd.DataFrame.from_dict(data_dict)
            df = df.append(comb_df)
            print(len(df))
        except Exception as e:
            print(e)
            print('%s deletions are lethal' % mut_name)
            continue
    return df

## MODIFY FOR READING INSIDE FRAMEWORK DIRECTORY
def analyse_results(set_up_params, model, target_reaction, sorting_param, analysis_type='gene_candidates'):
    #get the main features of results and compute carbon yield for every strategy found in OptKnock:
    reaction_set = set()
    count_dict = {rxn.id:0 for rxn in model.reactions}
    flux_dict = {rxn.id:[] for rxn in model.reactions}
    biomass_dict = {rxn.id:[] for rxn in model.reactions}

    strategy_list = []
    flux_list = []
    biomass_list = []
    c_list = []
    return_var = []

    for c_l in set_up_params['max_cl_range']:
        for m_ko in set_up_params['max_knock_out_range']:
            for r in range(set_up_params['replicates']):
                try:
                    experiment = '_'.join(set_up_params['project'].split('_')[0:-2])
                    filename= experiment+'_'+str(m_ko)+'_ko_'+str(c_l)+'_'+str(r)+'_rep.csv'
                    filepath= '/'.join([set_up_params['dir_path'], set_up_params['project'], set_up_params['framework'], filename])
                    print('Reading %s' % filepath)
                    df = pd.read_csv(filepath)

                    for strategy in df.reactions.values.tolist():
                        if analysis_type == 'strategies_eval' or 'strategies_eval' in analysis_type:
                            strategy_list.append(list(ast.literal_eval(strategy)))
                            flux_list.append(df.loc[df['reactions']== strategy, target_reaction].values[0])
                            biomass_list.append(df.loc[df['reactions']== strategy, 'biomass'].values[0])
                            c_list.append(c_l)

                        if analysis_type == 'gene_candidates' or 'gene_candidates' in analysis_type:
                            strategy_set = ast.literal_eval(strategy)
                            flux_yield = df.loc[df['reactions']== strategy, target_reaction].values[0]
                            biomass = df.loc[df['reactions']== strategy, 'biomass'].values[0]
                            reaction_set |= strategy_set
                            print(reaction_set)
                            for rxn in strategy_set:
                                count_dict[rxn] += 1
                                flux_dict[rxn].append(flux_yield)
                                biomass_dict[rxn].append(biomass)

                except Exception as e:
                    print(e)
                    print('%s has no data!' % filepath)
                    continue

    if analysis_type == 'strategies_eval' or 'strategies_eval' in analysis_type:
        raw_data = {'Strategy': strategy_list,
                    'Flux': flux_list,
                    'Biomass': biomass_list,
                    'N_of_deletions': [len(s) for s in strategy_list],
                    'C_limit': c_list}
        raw_df = pd.DataFrame.from_dict(raw_data)
        raw_df.to_csv('/'.join([set_up_params['dir_path'], set_up_params['project'], set_up_params['framework'], experiment+'_optknock_strategy_analysis.csv']), index=False)
        return_var.append(raw_df)

    if analysis_type == 'gene_candidates' or 'gene_candidates' in analysis_type:
        max_biomass = max([b for b_l in biomass_dict.values() for b in b_l])
        max_flux = max([f for f_l in flux_dict.values() for f in f_l])
        reactions_in_strategies = len({k:v for k,v in count_dict.items() if count_dict[k] > 0})
        final_candidate_data = {'Reaction' : [r 
                                              for r in count_dict.keys() if count_dict[r] > 0],
                                'Presence' : [count_dict[r]/reactions_in_strategies
                                              for r in count_dict.keys() if count_dict[r] > 0],
                                'Normalised_max_Flux' :[max(flux_dict[r])/max_flux
                                                         for r in flux_dict.keys() if len(flux_dict[r]) > 0],
                                'Normalised_max_biomass' : [max(biomass_dict[r])/max_biomass 
                                                             for r in biomass_dict.keys() if len(biomass_dict[r]) > 0]}


        results_df = pd.DataFrame.from_dict(final_candidate_data).sort_values(by=sorting_param, ascending=False)
        results_df.to_csv('/'.join([set_up_params['dir_path'], set_up_params['project'], set_up_params['framework'],experiment+'_optknock_gene_analysis.csv']), index=False)
        return_var.append(results_df)

    if len(return_var) == 1:
        return_var = return_var[0]

    return return_var

#To avoid the incorrect filtering of blocked reactions due to
#diverging fluxes towards non-objective biomass reactions, a
#function will be defined which purge all biomass but the target
#one from the model. For this, we are going to asume that biomas  
#reactions meet thefollowing criteria:
#    1. They have non integer coefficients for reaction reactants
#    2. They are the top reactions concerning the nº of metabolites
#After deleting all biomass reactions, the function will check
#if the model is feasible and if it is not, it will raise a 
#warning and return the given model
def purge_non_objective_biomass(model, target_biomass, n_of_biomass_reactions=None, threshold=12):
    new_model = model.copy()
    target_biomass_rxn = new_model.reactions.get_by_id(target_biomass)
    reactions_with_nonint_coeffs = [r for r in new_model.reactions 
                                    if all([not coeff.is_integer() for coeff in r.metabolites.values() if coeff<0])]
    
    
    n_of_metabolites_dict = { r.id : len(r.metabolites) for r in reactions_with_nonint_coeffs }
    n_of_metabolites_dict = dict(sorted(n_of_metabolites_dict.items(), key=lambda item: item[1], reverse=True))

    
    if n_of_biomass_reactions is None:
        #if no number of biomass reactions is given, filtering of reactions will be
        #done by a defined thershold
        n_of_biomass_reactions = threshold
        
    biomass_reactions = list(n_of_metabolites_dict.keys())[:n_of_biomass_reactions]
            
    for r in biomass_reactions:
        if r != target_biomass :
            #assure the reaction is not an intermediary biomass reaction supporting
            #the BOF with reactants
            rxn_products = [p.id for p in new_model.reactions.get_by_id(r).products]
            if not any([s.id in rxn_products for s in target_biomass_rxn.reactants]):
                new_model.reactions.get_by_id(r).bounds = (0,0)
            
    if new_model.optimize().status == 'optimal':
    
        return new_model
    
    else:
        print('WARNING: THE REACTION PURGING YIELDS AN INFEASIBLE MODEL, PLEASE PURGE THEM MANUALLY')
        print('RETURNING INPUT MODEL...')
        return model
    
#GENERATE USEFUL DATA FROM FSEOF RESULT
def process_fseof_result(fseof_result):
    data = fseof_result.data_frame
    data_points = len(data.columns)
    #1.ACCOUNT FOR FLUX DIFFERENCE
    data['flux_difference'] = [row[data_points]-row[1] for index, row in data.iterrows()]
    #2.REMOVE THOSE WITH NO CHANGE
    data = data.loc[data['flux_difference']!=0]
    #3.ORDER REACTION BY CHANGES IN FLUX
    data.sort_values(by='flux_difference', ascending=False)
    #4.GET ONLY DIRECTLY COUPLED REACTIONS:
    data['directed_coupling'] = [(all([ (abs(row[i+1])-abs(row[i])) > 0 for i in range(1,data_points)]) or 
                                  all([ (abs(row[i+1])-abs(row[i])) < 0 for i in range(1,data_points)]))
                                 for index, row in data.iterrows() ]
    
    data['regulation_type'] = ['U' if (row[1]<row[data_points] and row[1]>=0)
                               or (row[1]>row[data_points] and row[1]<0)
                               else 'D'
                               for index, row in data.iterrows()]

    data = data.sort_values(by='flux_difference', ascending=False)
    print('Saving proccessed result in "fseof_analysis.csv"...')
    data.to_csv('fseof_analysis.csv')
    return data

def plot_fseof_analysis(data, target_reaction, threshold,  regulation_type='U', to_remain=[]):
    
    target_row = pd.DataFrame([data.loc[target_reaction.id]])
    df = data.loc[(data['flux_difference']<=-threshold) | (data['flux_difference']>=threshold) ]
    #filter by the type of regulation you are interested in
    df = df.loc[df['regulation_type']==regulation_type]
    df= pd.concat([df, target_row])
    df.pop('flux_difference')
    df.pop('directed_coupling')
    df = df.T
    #to avoid duplication of column number in the case that target reaction roww meets the first criteria:
    if df.columns.tolist().count(target_reaction.id) > 1:
        df = df.loc[:,~df.columns.duplicated()]
    
    if len(to_remain) > 0:
        df = df[df.columns.intersection(to_remain)]
    display(df)
    fig = px.line(df, x=target_reaction.id, y=list(set(df.columns.tolist())-set(target_reaction.id)))

    fig.update_layout(    
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        title=dict(text='FSEOF result for '+''.join([m.name for m in target_reaction.metabolites])+' production', x=0.5),
        xaxis=dict(title=' '.join(['Flux',target_reaction.id]), showline=True, linewidth=3.0, linecolor='black',
                   ticks="inside", tickwidth=3.0, tickcolor='black', ticklen=5),
        yaxis=dict(title='Flux of top coupled reactions', showline=True, linewidth=3.0, linecolor='black',
                   ticks="inside", tickwidth=3.0, tickcolor='black', ticklen=5),
        font=dict(size=30),
        legend=dict(title='Reactions', font=dict(size=16))
    )

    fig.update_traces(marker=dict(line=dict(width=3.0),opacity=0.8))
    
    filename = 'FSEOF_'+target_reaction.id+'.svg'
    fig.show()
   
    print('Saving image result as %s ...' % filename)
    fig.write_image(filename)
    return df.columns.tolist()