% This script does all simulations 

addpath('simulation')
addpath('modeling')

%% IDENTIFIABILITY AND RECOVERY
%% 1. Simulate RL with broad parameters (used for identifiability and recovery)
%Simulate choice behavior for each model and do the "cross-fitting" 
% of all models (i.e. fit the data of each simulated model on all possible
% models). 

outfilename = 'SimsRLIdentRecov'; %filename for saving the simulation results
modelsinfoscript = 'modeling/loadModelsInfoAbsRelCon'; %script containing the specifications for the models 
outcomes = [0.1,1]; % outcome magnitudes to simulate  
iStartSim = 1;  %starting number of simulation, bigger than 1 if appending to existing simulations file
seed = randi(1000); %seed for pseudo-random number generation
whichmodels = []; %empty to automatically select all;
modelSimulateLearning(outfilename,'modeling/loadModelsInfoAbsRelCon',1,[0.1,1], iStartSim, seed,whichmodels)

%% 2. Simulate confidence based on RL simulation results (used for identifiabilty and recovery)
% Simulate and cross-fit confidence regressions within a given learning 
% model (the winning one), using the results (choice and values) simulated
% in the previous step. The script below also makes use of the file 
% containing the confidence regression coefficients for real participants,
% which it uses to define distributions from which to sample the generative
% confidence coefficients of the simulations (unlike for the RL parameters, 
% there were no clear broad priors informed by the literature).

modelSimulateConfidenceLogit;

%% 3. Run identifiability and parameter recovery
% Calculate BMC and correlations of generated and recovered RL parameters 
% and confidence coefficients, save results
modelIdentRecovChoice;
modelIdentRecovConfidence;

%% SIMULATIONS WITH DISTRIBUTIONS BASED ON PARTICIPANTS PARAMETERS
% These are used for testing the replicability of observed statistical 
% effects in real participants (Figures 3, 9 and S15, and S16).
modelSimulateSubjParamsConfLogit; %do simulations
simulationsNPReg; %perform non-parametric regressions of behavior over RL 
                  %parameters and confidence coefficients
