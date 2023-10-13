%%%Make matrices for model identifiability and parameter 
%recovery for confidence
%(Palminteri et al 2016,
%https://doi.org/10.1016/j.tics.2017.03.011)

clear
close all
clc
addpath('ModelingFuncs')

%%% get list of files. Each file has one simulation
% file = dir ('Results/SIMS_conf2022_04_08_CONFIDENCELOGIT*');
% load([file(1).folder,filesep,file(1).name])

load('Results/SimsConfIdentRecov.mat')

nsims = isim;
nsubs = size(regress.BIC,2);

f = load('Results\reg_conflogit_learning_dqabs.mat','whichLearnModel' ,'isConfPrev' ,'whichConfModel');
% whichLearnModel = f.whichLearnModel; isConfPrev = f.isConfPrev; whichConfModel = f.whichConfModel;
% idcmodels = 1:numel(whichLearnModel)%find(~modelvars.isConfPrev);
% whichmodelsLT = find((whichLearnModel == 11));
% whichmodelsTT = find((whichLearnModel == 11))
% whichmodelsLT = find((whichLearnModel == 11) &isConfPrev);
% whichmodelsTT = find((whichLearnModel == 11)  &~isConfPrev);

idcmodels = 1:16;
whichmodelsLT = 1:16;
whichmodelsLTReduced = 7:10;

whichmodelsTT = 1:16;
whichmodelsTTReduced = 1:4;

nmodelsLT = numel(whichmodelsLT);
nmodelsTT = numel(whichmodelsTT);
BICAllLT = nan(nsims,nsubs,nmodelsLT,nmodelsLT);
BICAllTT = nan(nsims,nsubs,nmodelsTT,nmodelsTT);

% prepare parameters for recovery
% find index and number of parameters for best model
% modelLT = find(whichLearnModel(idcmodels) == 11 & whichConfModel(idcmodels) == 2 & isConfPrev(idcmodels));
% modelTT = find(whichLearnModel(idcmodels) == 1 & whichConfModel(idcmodels) == 2 & ~isConfPrev(idcmodels));    
modelLT = 8;
modelTT = 2;
modelLT2 = 16;
modelTT2 = 14;

nparamsLT = numel(genparamsconf{1,1,modelLT})+1; %number of parameters of full model
nparamsTT = numel(genparamsconfpost{1,1,modelTT})+1; %number of parameters of full model

%%% get matrices generative and recovered parameters (sub x par x sim)
% genparLT = nan(nsubs,nparamsLT,nsims);
% recparLT = nan(nsubs,nparamsLT,nsims);
% genparTT = nan(nsubs,nparamsTT,nsims);
% recparTT = nan(nsubs,nparamsTT,nsims);
% 
% clear regress
% 
%% model identifiability LT
[pxpLT,bmLT] = doModelIdent(regress.BIC(1:isim,:,:,:),whichmodelsLT);
[pxpTT,bmTT] = doModelIdent(regress.BICPost(1:isim,:,:,:),whichmodelsTT);


[pxpLTReduced,bmLTReduced] = doModelIdent(regress.BIC(1:isim,:,:,:)./2,whichmodelsLTReduced);
[pxpTTReduced,bmTTReduced] = doModelIdent(regress.BICPost(1:isim,:,:,:)./2,whichmodelsTTReduced);

%% parameter recovery (for specific model)
for isub  = 1:nsubs
    for isim = 1:nsims
        genparLT(isub,:,isim) = genparamsconf{isim,isub,modelLT};
        recparLT(isub,:,isim) = regress.coeffs{isim,isub,modelLT,modelLT};
        genparTT(isub,:,isim) = genparamsconfpost{isim,isub,modelTT};
        recparTT(isub,:,isim) = regress.coeffsPost{isim,isub,modelTT,modelTT};
        
        genparLT2(isub,:,isim) = genparamsconf{isim,isub,modelLT2};
        recparLT2(isub,:,isim) = regress.coeffs{isim,isub,modelLT2,modelLT2};
        genparTT2(isub,:,isim) = genparamsconfpost{isim,isub,modelTT2};
        recparTT2(isub,:,isim) = regress.coeffsPost{isim,isub,modelTT2,modelTT2};
        
    end
end
    

[RestLT,R2estLT] = doParamRecov(genparLT,recparLT);
[RestTT,R2estTT] = doParamRecov(genparTT,recparTT);

[RestLT2,R2estLT2] = doParamRecov(genparLT2,recparLT2);
[RestTT2,R2estTT2] = doParamRecov(genparTT2,recparTT2);

infostr = struct('nsim',nsims,'whichmodelsLT',whichmodelsLT,'modelLT',modelLT,...
    'whichmodelsTT',whichmodelsTT,'modelTT',modelTT);

save('Results\ResultsIdentRecovConfidence.mat',...
    'bmLT','pxpLT','RestLT','R2estLT','genparLT','recparLT',...
    'bmTT','pxpTT','RestTT','R2estTT','genparTT','recparTT',...
    'pxpLTReduced','bmLTReduced','pxpTTReduced','bmTTReduced',...
    'genparLT2','recparLT2','genparTT2','recparTT2','RestLT2','R2estLT2',...
    'RestTT2','R2estTT2',...
    'infostr')

