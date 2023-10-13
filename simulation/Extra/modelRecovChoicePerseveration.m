%%%Make matrices for parameter
%recovery (Palminteri et al 2016,
%https://doi.org/10.1016/j.tics.2017.03.011)

clear
close all
clc
addpath('ModelingFuncs')
% load simulation data

% simfiles = {'Results\SIMS_conf2020_05_01.mat','Results\SIMS_50c.mat'};
% outfilenames = {'SIMS_ResultsIdentRecov.mat','SIMS_ResultsIdentRecov50c.mat'};

simfiles = {'Results\SIMS_ExtraModels'};
outfilenames = {'SIMS_ResultsRecovPerseveration.mat'};

for isimfile = 1
    fname = simfiles{isimfile};

    load(fname);
    if exist('isim','var') && isim~=nsims; %isim stores ongoing simulation if not finished
        nsims = isim-1;
    else
        nsims = nsims; %n_sims stores total number of simulations
    end
    nsubs = size(LAME,2);
    %% parameter recovery (for specific model)
    model = 1;
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
    
    infostr = struct('nsim',nsims,'nmodels',1,'model',model);
    
    % save('Results\SIMS_ResultsIdentRecov.mat','bm','pxp','Rest','R2est',...
    %     'genpar','recpar','infostr')
    
    save(['Results\',outfilenames{isimfile}],'Rest','R2est',...
        'genpar','recpar','infostr')
end