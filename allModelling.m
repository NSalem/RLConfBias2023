% This scripts does all the fitting of learning models and confidence
% regressions

addpath('modeling')
%% 1. Reinforcement Learning models
%fit learning models, saving parameters, model evidence, and 
% estimated learnt values for each participant
fitModels; 

%% 2. Confidence regressions
%Regress Learning task confidence on estimated learnt values for each
%participant, save results
regressConfTestLogit; 

%Regress Transfer task confidence on estimated learnt values for each
%participant, save results
regressConfPostLogit;

%Regress confidence from both tasks on estimated learnt values, using
%regression models with and without different parameters for Learning and Transfer
regressConfDual;
