%%% fit RELASYM perseveration model on simulations produced without
%%% perseveration

%load sims

clear
close all
clc
addpath('ModelingFuncs')

load('Results\SIMS_conf2022_11_11.mat')

if exist('isim','var') && isim~=nsims; %isim stores ongoing simulation if not finished
    nsims = isim-1;
else
    nsims = nsims; %n_sims stores total number of simulations
end
nsubs = size(LAME,2);
%% parameter recovery (for specific model)
modelgen = 1;
modelrec = 2;
nparams = numel(genparams{1,1,modelgen}); %number of parameters of full model
%get matrices generative and recovered parameters (sub x par x sim)
genpar = nan(nsubs,nparams,nsims);
recpar = nan(nsubs,nparams,nsims);

for isub  = 1:nsubs
    for isim = 1:nsims
        genpar(isub,:,isim) = genparams{isim,isub,modelgen};
        recpar(isub,:,isim) = parametersLPP{isim,isub,modelgen,modelrec}(:,1:nparams);
    end
end

[Rest,R2est] = doParamRecov(genpar,recpar);

infostr = struct('nsim',nsims,'nmodels',1,'model',2);

% save('Results\SIMS_ResultsIdentRecov.mat','bm','pxp','Rest','R2est',...
%     'genpar','recpar','infostr')

save('Results\SIMS_ResultsCrossRecovPerseveration','Rest','R2est',...
    'genpar','recpar','infostr')
