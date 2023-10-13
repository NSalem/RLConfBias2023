% Fit confidence models on hidden variables of leraning model for the
% transfer task
% Full confidence models: y~S(deltaQ+y(t-1)+bias), where S is sigmoid
% function, and 'bias' can be 0, Qc, sigmaQ, or V (the latter only works 
% if the RL model is contextual)


clear

addpath('ModelingFuncs\');
addpath('helperfuncs');
resultsdir = ['Results',filesep];
datadir = ['data',filesep];

learnmodels = 11;

%%% iterate for different implementations of difficulty variable
outfilenames = {'reg_conflogit_posttest','reg_conflogit_posttest_dqabs','reg_conflogit_posttest_nodq','reg_conflogit_posttest_dqgoodbad','reg_conflogit_posttest_pc'};
RLvarsfiles = {'RLVars','RLVars','RLVars','RLVars','RLVars'};

confmodels ={{'dQ','dQ+Qc','dQ+QcplusQu','dQ+V', 'dQ+Qc+V','dQ+QcplusQu+V'},... 
             {'dQabs','dQabs+Qc','dQabs+QcplusQu','dQabs+V', 'dQabs+Qc+V','dQabs+QcplusQu+V'},...
             {'Qc','QcplusQu','V','Qc+Qu','Qc+Qu+V'},... 
             {'dQGoodBad','dQGoodBad+Qc','dQGoodBad+QcplusQu','dQGoodBad+V', 'dQGoodBad+Qc+V','dQGoodBad+QcplusQu+V'},...
             {'pc','pc+Qc','pc+QcplusQu','pc+V', 'pc+Qc+V','pc+QcplusQu+V'},...
             };
%repeat with previous confidence 
for ifile = 1:numel(confmodels)
    confmodels{ifile} =[{confmodels{ifile}{:}},strcat({confmodels{ifile}{:}},'+ confprev')]; 
end

% fit and save 
for ifile = 1:numel(outfilenames)
    regressConfTT(['Results/',RLvarsfiles{ifile}],learnmodels,confmodels{ifile},['Results/',outfilenames{ifile}]);
end
