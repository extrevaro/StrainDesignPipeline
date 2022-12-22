load Irene2.mat
model = changeRxnBounds(model,{'EX_glc___D_e'},0,'l')
model = changeRxnBounds(model,{'EX_ac_e'},-10,'l')
model=changeObjective(model,'EX_btd_RR(e)')
sol = optimizeCbModel(model,'max','one')