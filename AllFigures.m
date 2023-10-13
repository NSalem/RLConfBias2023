%%%% Plot all figures used in paper %%%%%%

addpath('makeplots')
addpath('makeplots/AdditionalContextualization/')

%% Main results 

%%% Model fits
Figure_ModelFitsLT %Learning
Figure_ModelFitsTT %Transfer 
close all
%%% Model comparison choice
Figure_BMCchoice 
% Figure_BMCchoiceExtended (include V and V+Q context models)
close all

%%% Model comparison confidence
Figure_BMCconfLT
Figure_BMCconfTT
close all
% hidden vars LT
Figure_HiddenVarConfLT %plot hidden variables for learning task
% hidden vars TT
Figure_HiddenVarsTT
close all force
%%% Param distributions
Figure_ParamDistribs
close all force

%%%  Correlate fitted parameters with behavioral measures 
% (accuracy, confidence, overconfidence, valence bias).
% Requiers having ran simulationsNPReg.m (see allSimulations.m)
corrParamsBehav;

% Correlation of learning and confidnece parameters and confidence biases
% + multiple linear regression for valence bias and overconfidence on
% leraning params.
% Requires having ran corrParamsBehav.m

Figure_CorrParamsBehavIndivPlot;
close all force

% Model identifiability and parameter recovery
Figure_IdentRecovChoice
Figure_IdentRecovConfidence %recovery for best model AND Qc+Qu+V. It seems Qc and V are too correlated in transfer?
close all force
% Qc+Qu / Qc+Qu+V conf model param dist
Figure_ParamsQcQuV 
close all force

% Correlation of behavioral biases between tasks
Figure_CorrBehavLearnTransfer

% Correlation of confidence coefficients between tasks
Figure_CorrConfParamsLearnPost

%Bayesian Model Comparison of same vs different parameters in learning and
%transfer
Figure_BMC_ConfDual

close all force

%%% Stats 
Figure_StatsBehavior
Figure_StatsBehaviorTT 
close all force
%%% Difference between confidence models in transfer task (simulations)
Figure_ConfmodelsDiffTT;

%%% comparison prevconf vs no-prev-conf
Figure_ModelComparisonConfidencePrevVsNot

%%% comparison different difficulty models conf
Figure_ModelComparisonsConfidenceDifficulty 
close all force


%% Other
%%% legend of symbols (valence, info, pairing, utility...)
addpath('makeplots/Extra')
Figure_ConditionLegend
close all force
