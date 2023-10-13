%%%Compare confidence (which are nested within the learning models).
%%% Learning models: ABS, REL0, RELQu, RELother, RELlast
%%% Confidence models: y~deltaQ+y(t-1)+context
%%% where 'context' can be 0, Qc, sigmaQ, or V (the latter only if the RL model is contextual)

%%% For ABS and RELother models:
%%% 1- load hidden variables (Q, V) for each subject based on previously fitted parameters
%%% 2- For each context variant, regress confidence on the chosen set of variables
%%%(y~deltaQ+y(t-1)+context)

clear

addpath('ModelingFuncs\');
addpath('helperfuncs');
resultsdir = ['Results',filesep];
datadir = ['data',filesep];

rlnames = {'ROther2','ROtherNorm','Perseveration'};


for ilearnmodel = 1:3
    %     learnmodels = ilearnmodel;
    
    % outfilenames = {'reg_conflogit_learning','reg_conflogit_learning_splitbeta','reg_conflogit_learning_dqabs','reg_conflogit_learning_splitbeta_dqabs','reg_conflogit_learning_nodq','reg_conflogit_learning_splitbeta_nodq'};
    % RLvarsfiles = {'RLVars','RLVars_BetaSplit','RLVars','RLVars_BetaSplit','RLVars','RLVars_BetaSplit'};
    % difficultyTerms = {'+dQ','+dQ','+dQabs','+dQabs','',''};
    
    
    outfilenames = {['reg_conflogit',rlnames{ilearnmodel},'_learning_dqabs']};
    RLvarsfiles = {'RLvars_ExtraContext'};
    
    confmodels ={{'dQabs','dQabs+Qc','dQabs+QcplusQu','dQabs+V', 'dQabs+Qc+V','dQabs+QcplusQu+V'}};
    
    %repeat with previous confidence
    for ifile = 1
        confmodels{ifile} =[{confmodels{ifile}{:}},strcat({confmodels{ifile}{:}},'+ confprev')];
    end
    
    % fit and save
    for ifile = 1
        regressConfLT(['Results/',RLvarsfiles{ifile}],ilearnmodel,confmodels{ifile},['Results/',outfilenames{ifile}]);
    end
    
end