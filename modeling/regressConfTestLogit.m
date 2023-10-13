% Fit confidence models (regression) on hidden variables of leraning 
% model for the learning task.
% Full confidence models: y~S(deltaQ+y(t-1)+bias), where S is sigmoid
% function, and 'bias' can be 0, Qc, sigmaQ, or V (the latter only works if 
% the RL model is contextual)


clear
addpath('ModelingFuncs\');
addpath('helperfuncs');
resultsdir = ['Results',filesep];
datadir = ['data',filesep];

learnmodels = 11; %11: ASYMRELOther

% contextTerms = {'','+ Qc','+ QcplusQuc','+ V','+ Qc + Quc','+ Qc + Quc + V'}; %

%%% iterate for different "difficulty/discriminability" variable

outfilenames = {'reg_conflogit_learning_dqabs','reg_conflogit_learning','reg_conflogit_learning_nodq','reg_conflogit_learning_dqgoodbad','reg_conflogit_learning_pc'};
RLvarsfiles = {'RLVars','RLVars','RLVars','RLVars','RLVars'};

confmodels ={{'dQabs','dQabs+Qc','dQabs+QcplusQu','dQabs+V', 'dQabs+Qc+V','dQabs+QcplusQu+V'},...
            {'dQ','dQ+Qc','dQ+QcplusQu','dQ+V', 'dQ+Qc+V','dQ+QcplusQu+V'},... 
             {'Qc','QcplusQu','V','Qc+Qu','Qc+Qu+V'},... 
             {'dQGoodBad','dQGoodBad+Qc','dQGoodBad+QcplusQu','dQGoodBad+V', 'dQGoodBad+Qc+V','dQGoodBad+QcplusQu+V'},...
             {'pc','pc+Qc','pc+QcplusQu','pc+V', 'pc+Qc+V','pc+QcplusQu+V'},...
             };

%repeat with previous confidence 
for ifile = 1:numel(confmodels)
    confmodels{ifile} =[{confmodels{ifile}{:}},strcat({confmodels{ifile}{:}},'+ confprev')]; 
end

% fit and save 
for ifile = 1:numel(confmodels)
    regressConfLT(['Results/',RLvarsfiles{ifile}],learnmodels,confmodels{ifile},['Results/',outfilenames{ifile}]);
end
