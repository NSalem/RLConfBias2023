%% fit models and save parameters %%

addpath('ModelingFuncs\');
addpath('helperfuncs\');
addpath('modeling')

%%% fit and save results for all 10 exps
rng(0)
modelFitPARALLEL('Results\data_all', 'CON2','loadModelsInfoAbsRelCon.m',[1:10],[],[], 1,1);
saveModelVars('Results\data_all','Results\RLMODELCON2exp','RLvars_all')

%%%% fit and save results for confidence exps only
rng(0)
modelFitPARALLEL('Results\data_all_CONF', 'CON','loadModelsInfoAbsRelCon.m',[],[],[], 1,1);
saveModelVars('Results\data_all_CONF','Results\RLMODELCONexp','RLvars')

%%% fit extra models
% rng(0)
% modelFitPARALLEL('Results\data_all', 'CON2_ExtraContext','modeling/Extra/loadModelsInfoAdditionalContextualization.m',[],[],[], 1,1);
% saveModelVars('Results\data_all','Results\RLMODELCON2_ExtraContextexp','RLvars_ExtraContext_all')

%%% fit extra models for confidence exps only
% rng(0)
% modelFitPARALLEL('Results\data_all_CONF', 'CON_ExtraContext','loadModelsInfoAdditionalContextualization.m',[],[],[], 1,1);
% saveModelVars('Results\data_all_CONF','Results\RLMODELCON_ExtraContextexp','RLvars_ExtraContext')
