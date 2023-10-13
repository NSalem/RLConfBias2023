%%%Make matrices for model identifiability and parameter
%recovery (Palminteri et al 2016,
%https://doi.org/10.1016/j.tics.2017.03.011)

clear
close all
clc
addpath('ModelingFuncs')
% load simulation data

% simfiles = {'Results\SIMS_conf2020_05_01.mat','Results\SIMS_50c.mat'};
% outfilenames = {'SIMS_ResultsIdentRecov.mat','SIMS_ResultsIdentRecov50c.mat'};

simfiles = {'Results\SimsRLIdentRecov'};
outfilenames = {'ResultsIdentRecovChoice'};


for isimfile = 1
    fname = simfiles{isimfile};
    % fname = 'Results\SIMS_conf2020_05_01.mat';
    % fname = 'Results\SIMS_reducedTrials.mat';
    
    load(fname);
    
    if exist('isim','var') && isim~=nsims; %isim stores ongoing simulation if not finished
        nsims = isim-1;
    end
    
    nsubs = size(LAME,2);
    whichmodels = 1:14;
    nmodels = numel(whichmodels);
    
    %% model identifiability
    criterion = LAME;
    [pxp,bm] = doModelIdent(LAME(1:nsims,:,:,:),whichmodels);
    
    [pxpReduced,bmReduced] = doModelIdent(LAME(1:nsims,:,:,:),[1:5,8:12]);

    %% parameter recovery (for specific model)
    model = 11;
    nparams = numel(genparams{1,1,model}); %number of parameters of full model
    %get matrices generative and recovered parameters (sub x par x sim)
    genpar = nan(nsubs,nparams,nsims);
    recpar = nan(nsubs,nparams,nsims);
    
    for isub  = 1:nsubs
        for isim = 1:nsims
            genpar(isub,:,isim) = genparams{isim,isub,model};
            recpar(isub,:,isim) = parametersLPP{isim,isub,model,model};
        end
    end
    
    [Rest,R2est] = doParamRecov(genpar,recpar);
    
    infostr = struct('nsim',nsims,'nmodels',nmodels,'model',model);
    
    % save('Results\SIMS_ResultsIdentRecov.mat','bm','pxp','Rest','R2est',...
    %     'genpar','recpar','infostr')
    
    save(['Results\',outfilenames{isimfile}],'bm','pxp','pxpReduced','bmReduced','Rest','R2est',...
        'genpar','recpar','infostr')
end