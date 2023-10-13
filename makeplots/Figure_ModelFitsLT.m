% Plot behavior and model fits (based
% on best learning+confidence models) for Learning Task. 
% Plots incluided: Learning and confidence curves, violinplots of
% participant averages.

clear;

%% Load data
% load('Results\data_all.mat');
addpath('helperfuncs');
rlVars = load('Results\RLvars_all.mat');
imodel = 11;
reg = load('Results\reg_conflogit_learning_dqabs.mat');
pcMod = squeeze(nanmean(rlVars.pc(imodel,:,1:4,:,:),3)).*100;

idc = find(reg.isConfPrev & reg.whichLearnModel==imodel & strcmp(reg.confBias, '+Qc'));
confMod = squeeze(nanmean(reg.confmodel(:,idc,:,:,:),3));
clear reg


%% Average session data
corr_mat = squeeze(nanmean(rlVars.correct,2));
conf_mat =  squeeze(nanmean(rlVars.conf,2));
ntrials = size(corr_mat,3);

%% Plot results
%%% acc for all data
[h1,h2] = makePlotsModelFitsLT(corr_mat, [], pcMod,[]);
% saveas(h1,'Plots/modelFitsLearnAccuracyAllData_Time.svg');
% saveas(h2,'Plots/modelFitsLearnAccuracyAllData_Violins.svg');

%%% confidence exps only (acc, conf, overconf)%%%
[h1,h2] = makePlotsModelFitsLT(corr_mat, conf_mat, pcMod,confMod);
% saveas(h1,'Plots/modelFitsLearnAccuracyAllData_Time.svg');
% saveas(h2,'Plots/modelFitsLearnAccuracyAllData_Violins.svg');


