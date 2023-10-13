function modelSimulateLearning(outfilename,modelsinfoscript,useParallel,outcomes, iStartSim, seed,whichmodels)
% Simulate choice (learning) behavior and hidden variables (values), and
% get fits across models. 
% Populations are generated for each model, and their resulting 
% behavior (choices) is fitted on each model. Results are saved to a file
% Parameters are drawn from the same distributions used for priors in
% fitting. 

%INPUT
%outfilename: filename to save results
%modelsinfoscript: path to script defining specifications of different
%models (if not specified or empty a UI dialog opens to select one)
%useParallel: bool, whether to use parallel computing for the fitting
%outcomes: size of outcomes used (vector of 2 elements), default = [0.1,1]
%iStartSim: number of starting simulation (default = 1, another number can
%be specified to add simulations to existing data file)
%seed: random seed to be used
%whichmodels: vector of which models to use from modelsinfoscript (default
%= all)

if nargin<1 || isempty(outfilename)
    outfilename = ['SimsRLIdentRecov_',datestr(date,'yyyy_mm_dd')];
end
if nargin<2 || isempty (modelsinfoscript)
    [file,path] = uigetfile;
    modelsinfoscript = [path,file];
end
addpath('ModelingFuncs\');
addpath ('helperfuncs\');

if nargin<3 || isempty(useParallel)
    useParallel = true;
end

if nargin <4 || isempty(outcomes)
    outcomes =[0.1,1];
end
if nargin<5 || isempty(outcomes)
    iStartSim = 1;
end

if nargin<6 || isempty(seed)
    seed = randi(1000);
end

if nargin <7
    whichmodels = [];
end

rng(seed);

nsims      = 20;
nsub       = 90;
nsess      = 3;
ntrials = 24;
nreppost = 4;
ntrialspost = nchoosek(8,2)*nreppost;
% usepost = 1;

if useParallel
    if isempty(gcp('nocreate'))
        parpool
    end
    parArg = Inf;
else
    parArg = 0;
end

%% Functions of distribution of generative RL parameters
gen.beta  = @()random('Gamma',1.2,5);
gen.betaTT  = @()random('Gamma',1.2,5);
gen.lr1 = @()random('Beta',1.1,1.1);
gen.lr2 =@() random('Beta',1.1,1.1);
gen.lr3 = @()random('Beta',1.1,1.1);
gen.ROtherWeight = @()random('Beta',1.1,1.1);
gen.tau = gen.lr1;
gen.fi = @() random('Normal',0,1);

%% load models specifications
run(modelsinfoscript);

if ~isempty(whichmodels)
    modelsinfo = modelsinfo(whichmodels);
end

clear priorfuncs

%% initialize arrays, etc
nlearnmodels = (numel(modelsinfo));
if iStartSim ==1 %if starting new simulations
    genparams = cell(nsims,nsub,nlearnmodels);
    parametersLPP = cell(nsims,nsub,nlearnmodels,nlearnmodels);
    LAME = nan(nsims,nsub,nlearnmodels,nlearnmodels);
    LPP = nan(nsims, nsub, nlearnmodels,nlearnmodels); %LPP matrix of generative models by recovered models for each subj and sim
    ll = nan(nsims, nsub, nlearnmodels,nlearnmodels);
    %     aic= nan(n_sims, n_sub, nlearnmodels,nlearnmodels);
    %     bic = nan(n_sims, n_sub, nlearnmodels,nlearnmodels);
    choice_all = nan(nsims,nsub,nlearnmodels,nsess,ntrials*4);
    pc_all = nan(nsims,nsub,nlearnmodels,nsess,ntrials*4);
    s_all=nan(nsims,nsub,nlearnmodels,nsess,ntrials*4);
    r_all=nan(nsims,nsub,nlearnmodels,nsess,ntrials*4);
    c_all=nan(nsims,nsub,nlearnmodels,nsess,ntrials*4);
    choicePost_all=nan(nsims,nsub,nlearnmodels,ntrialspost);
    ss_all = nan(nsims,nsub,nlearnmodels,2,ntrialspost);
    
else %if continuing on an existing file
    load(['Results',filesep,outfilename,'.mat']);
end

tic
dummodelsinfo = modelsinfo;

load('Results\OldSims\SIMS_ResultsIdentRecov_2022_02_28.mat','genpar')
OUT = Generate_Outcomes(24,outcomes); % generate outcomes
for isim = iStartSim:nsims
    disp(['Simulation ', num2str(isim)]);
    %%%genearte parameters
    genparams_this_sim = getlearnparams(dummodelsinfo,nsub,gen);
    genparams(isim,:,:) = genparams_this_sim;
    %%%simulate action and hidden variables for each subject and learning
    %%%model
    [outSim] = simulateRLPop(genparams_this_sim,modelsinfo,nsess,ntrials,outcomes);
    
    %     [pc_all_sim,choice_all_sim,choicePost_all_sim,s_all_sim,ss_all_sim,r_all_sim,c_all_sim]=simulateRLPop(genparams_this_sim,dummodelsinfo,n_sess,ntrials,outcomes,nreppost);
    
    %%%fit models to generated actions
    [LAME_sim,parametersLPP_sim,ll_sim,LPP_sim] = fitModels(dummodelsinfo,nsub,nsess,ntrials,...
        outSim.s,outSim.ss,outSim.choice,outSim.choicePost,outSim.r,outSim.c);
    
    %%% MAKE MATRICES FOR ALL SIMS %%%
    choice_all(isim,:,:,:,:) = outSim.choice;
    pc_all(isim,:,:,:,:) = outSim.pc;
    s_all(isim,:,:,:,:) = outSim.s;
    r_all(isim,:,:,:,:) = outSim.r;
    c_all(isim,:,:,:,:) = outSim.c;
    choicePost_all(isim,:,:,:) = outSim.choicePost;
    ss_all(isim,:,:,:,:) = outSim.ss;
    LAME(isim,:,:,:) = LAME_sim;
    parametersLPP(isim,:,:,:) = parametersLPP_sim;
    LPP(isim,:,:,:) = LPP_sim;
    
    if ~exist(['Results',filesep,outfilename,'.mat'],'file')
        save(['Results',filesep,outfilename,'.mat'], 'LPP', 'LAME','parametersLPP',...
            'genparams','pc_all','choice_all','choicePost_all',...
            's_all','r_all','c_all','ss_all','ntrials', 'nsims', 'modelsinfo','isim')
    else
        save(['Results',filesep,outfilename,'.mat'], 'LPP', 'LAME','parametersLPP',...
            'genparams','pc_all','choice_all','choicePost_all',...
            's_all','r_all','c_all','ss_all','ntrials', 'nsims', 'modelsinfo','isim','-append')
    end
    
end
runtime = toc;
save(['Results',filesep,outfilename,'.mat'], 'LPP', 'LAME','parametersLPP',...
    'genparams','pc_all','choice_all','choicePost_all',...
    's_all','r_all','c_all','ss_all','ntrials', 'nsims', 'modelsinfo','isim','runtime','seed')

end


function [LAME,parametersLPP,ll,LPP] = fitModels(modelsinfo,n_sub,n_sess,n_trials,s_all, ss_all,a_all,aa_all,r_all,c_all)
nlearnmodels = numel(modelsinfo);
parametersLPP = cell(n_sub,nlearnmodels,nlearnmodels);
LAME = nan(n_sub,nlearnmodels,nlearnmodels);
LPP = nan(n_sub, nlearnmodels,nlearnmodels); %LPP matrix of generative models by recovered models for each subj and sim
ll = nan(n_sub, nlearnmodels,nlearnmodels);
% aic= nan(n_sub, nlearnmodels,nlearnmodels);
% bic = nan(n_sub, nlearnmodels,nlearnmodels);

for igenlearnmodel = 1:nlearnmodels
    disp(['Gen Model ',num2str(igenlearnmodel)])
    for isub = 1:n_sub
        s = squeeze(s_all(isub,igenlearnmodel,:,:));
        a = squeeze(a_all(isub,igenlearnmodel,:,:));
        ss = squeeze(ss_all(isub,igenlearnmodel,:,:))';
        aa = squeeze(aa_all(isub,igenlearnmodel,:));
        r = squeeze(r_all(isub,igenlearnmodel,:,:));
        c  = squeeze(c_all(isub,igenlearnmodel,:,:));
        
        ncond = numel(unique(s(:)));
        con_sub = nan(numel(s),1);
        cho_sub =  nan(numel(s),1);
        out_sub =  nan(numel(s),1);
        cou_sub =  nan(numel(s),1);
        for isess  = 1:n_sess
            k_in = (isess-1)*n_trials*ncond+1; %%% ncond not deffined!!!
            k_out =  isess*n_trials*ncond;
            con_sub(k_in:k_out,1) = (isess-1)*4 + s(isess,:);
            cho_sub(k_in:k_out,1) = a(isess,:);
            out_sub(k_in:k_out,1) = r(isess,:);
            cou_sub(k_in:k_out,1) = c(isess,:);
        end
        
        parfor ireclearnmodel =1:numel(modelsinfo)
            
            %%% do fits, get params for RL
            lb = modelsinfo{ireclearnmodel}.lb;
            ub = modelsinfo{ireclearnmodel}.ub;
            x0 = modelsinfo{ireclearnmodel}.x0;
            options = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 10000); % These increase the number of iterations to ensure the convergence
            
            %% fit recovered model %%%%
            % LPP (Laplace appriximation of the posterior probability) optimization
            [theseParamsLPP,thisLPP,~,~,~,~,hessian]=fmincon(@(x) GetModelLL_QLearner(x,modelsinfo{ireclearnmodel},con_sub,cho_sub,out_sub,cou_sub,aa,ss,1),x0,[],[],[],[],lb,ub,[],options);
            
            parametersLPP{isub,igenlearnmodel,ireclearnmodel} = theseParamsLPP;
            
            %                 parametersLPP{isim,isub,igenlearnmodel,ireclearnmodel} = theseParamsLPP;
            LPP(isub,igenlearnmodel,ireclearnmodel) = thisLPP;
            k = numel(modelsinfo{ireclearnmodel}.paramnames);
            LAME(isub,igenlearnmodel,ireclearnmodel) =  thisLPP - k/2*log(2*pi) + real(log(det(hessian))/2);%Laplace-approximated model evidence
            this_ll = GetModelLL_QLearner(theseParamsLPP,modelsinfo{ireclearnmodel},con_sub,cho_sub,out_sub,cou_sub,aa,ss,0);
            %             ll_test = GetModelLL_QLearner(theseParamsLPP,modelsinfo{ireclearnmodel},con_sub,cho_sub,out_sub,cou_sub,[],[],0);
            %             ll_post = this_ll-ll_test;
            ll(isub,igenlearnmodel,ireclearnmodel) = this_ll;
            %             bic(isub,igenlearnmodel,ireclearnmodel)=-2*-this_ll+k*log(n_sess*ntrials*4 + numel(aa));
            %             aic(isub,igenlearnmodel,ireclearnmodel)=-2*-this_ll+k; % l2 is already positive
        end
    end
end
end
