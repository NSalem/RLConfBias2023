% Simulate confidence behavior using previously simulated values (see
% modelSimulateLearning.m) and get fits across models. 
% Populations are generated for each model, and the resulting 
% behavior (confidence) is fitted on each model. Results are saved to a file
% Parameters are drawn from distributions fitted on our human sample.  

addpath('ModelingFuncs\')
iStartSim = 1;
doPostGen = 1;
%% Load RL sims file

RLsimfile = 'Results\SimsRLIdentRecov.mat';
outfilename = 'SimsConfIdentRecov';
load(RLsimfile,'choice_all','choicePost_all','s_all','ss_all','r_all','c_all','modelsinfo','parametersLPP','genparams')

whichlearnmodel = 11;

nsims = size(parametersLPP,1);
nsub = size(parametersLPP,2);
nsess = 3;
ntrials = 24;
ntrialspost = 112;

%% specify confidence models and generative parameter distributions
%
f = load('Results\reg_conflogit_learning_dqabs.mat','confcoeff','confRMSE','confBias','whichLearnModel' ,'isConfPrev' ,'whichConfModel','glmCONF','formulas');
f2 = load('Results\reg_conflogit_learning_nodq.mat','confcoeff','confRMSE','confBias','whichLearnModel' ,'isConfPrev' ,'whichConfModel','glmCONF','formulas');
whichLearnModel = [f.whichLearnModel,f2.whichLearnModel]; isConfPrev = [f.isConfPrev,f2.isConfPrev];
confBias = [f.confBias,f2.confBias([4,5,9,10])];
confcoeffLT = [f.confcoeff;f2.confcoeff([4,5,9,10])]; confRMSE = [f.confRMSE;f2.confRMSE([4,5,9,10])];
glmAll = [f.glmCONF(1,:),f2.glmCONF(1,[4,5,9,10])];

idcmodels = 1:size(confcoeffLT,1);%use all models from regression
modelvars = struct('whichLearnModel',whichLearnModel(idcmodels),'confBias',confBias(idcmodels),'isConfPrev',isConfPrev(idcmodels));
nmodels = numel(idcmodels);

coeffNames =cell(numel(idcmodels),1);
for iconfmodel = 1:numel(idcmodels)
    coeffNames{iconfmodel} = glmAll{1,iconfmodel}.CoefficientNames;
end
clear glmAll;

regpost = load('Results\reg_conflogit_posttest_dqabs.mat','confcoeff','confRMSE','whichLearnModel' ,'isConfPrev' ,'confBias','glmCONF','formulas');
regpost2 = load('Results\reg_conflogit_posttest_nodq.mat','confcoeff','confRMSE','whichLearnModel' ,'isConfPrev' ,'confBias','glmCONF','formulas');

confcoeffTT = [regpost.confcoeff';regpost2.confcoeff([4,5,9,10])']; confRMSEPost = [regpost.confRMSE';regpost2.confRMSE([4,5,9,10])'];


confmodelsinfo = cell(1,1);

confmodelsinfo{1} = struct('varnames',{{'intercept','dQabs'}},'ilearnmodel',1);
confmodelsinfo{2} = struct('varnames',{{'intercept','dQabs','Qc'}},'ilearnmodel',1);
confmodelsinfo{3} = struct('varnames',{{'intercept','dQabs','sigmaQ'}},'ilearnmodel',1);
confmodelsinfo{4} = struct('varnames',{{'intercept','dQabs','V'}},'ilearnmodel',1);
confmodelsinfo{5} = struct('varnames',{{'intercept','dQabs','Qc','V'}},'ilearnmodel',1);
confmodelsinfo{6} = struct('varnames',{{'intercept','dQabs','sigmaQ','V'}},'ilearnmodel',1);

confmodelsinfo{7} = struct('varnames',{{'intercept','dQabs','confPrev'}},'ilearnmodel',1);
confmodelsinfo{8} = struct('varnames',{{'intercept','dQabs','Qc','confPrev'}},'ilearnmodel',1);
confmodelsinfo{9} = struct('varnames',{{'intercept','dQabs','sigmaQ','confPrev'}},'ilearnmodel',1);
confmodelsinfo{10} = struct('varnames',{{'intercept','dQabs','V','confPrev'}},'ilearnmodel',1);
confmodelsinfo{11} = struct('varnames',{{'intercept','dQabs','Qc','V','confPrev'}},'ilearnmodel',1);
confmodelsinfo{12} = struct('varnames',{{'intercept','dQabs','sigmaQ','V','confPrev'}},'ilearnmodel',1);

confmodelsinfo{13} = struct('varnames', {{'intercept','Qc','Qu'}},'ilearnmodel',1);
confmodelsinfo{14} = struct('varnames', {{'intercept','Qc','Qu','V'}},'ilearnmodel',1);
confmodelsinfo{15} = struct('varnames', {{'intercept','Qc','Qu','confPrev'}},'ilearnmodel',1);
confmodelsinfo{16} = struct('varnames', {{'intercept','Qc','Qu','V','confPrev'}},'ilearnmodel',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


genconfdistsLT = cell(numel(confmodelsinfo),1);
genconfdistsTT = cell(numel(confmodelsinfo),1);

nmodels = numel(confmodelsinfo);

if iStartSim == 1
    genparamsconf = cell(nsims,nsub,nmodels);
    genparamsconfpost = cell(nsims,nsub,nmodels);
    gennoiseconf = cell(nsims,nsub,nmodels);
    gennoiseconfpost = cell(nsims,nsub,nmodels);
    regress.LL = nan(nsims, nsub, nmodels,nmodels);
    regress.AIC = nan(nsims, nsub, nmodels,nmodels);
    regress.BIC = nan(nsims, nsub, nmodels,nmodels);
    regress.coeffs = cell(nsims, nsub, nmodels,nmodels);
    regress.RMSE = nan(nsims, nsub, nmodels,nmodels);
    %     regress.BICCombTest = nan(nsims, nsub, nmodels,nmodels);
    %     regress.BICCombPost  = nan(nsims, nsub, nmodels,nmodels);
    %     regress.BICCombAll  = nan(nsims, nsub, nmodels,nmodels);
    %     regress.LLPost = nan(nsims, nsub, nmodels,nmodels);
    %     regress.AICPost = nan(nsims, nsub, nmodels,nmodels);
    regress.BICPost = nan(nsims, nsub, nmodels,nmodels);
    regress.coeffsPost = cell(nsims, nsub, nmodels,nmodels);
    regress.RMSEPost = zeros(nsims, nsub, nmodels,nmodels);
    %     coeffs  = cell(nsims, nsub, nmodels,nmodels);
    %     coeffsPost = cell(nsims, nsub, nmodels,nmodels);
    conf_all = nan(nsims,nsub,nmodels,nsess,ntrials*4);
    confPost_all = nan(nsims,nsub,nmodels,ntrialspost);
else
    load(['Results/',outfilename,'.mat'])
end


%% generate conf params
for isim = 1:nsims
    
    for igenmodel = 1:numel(confmodelsinfo)
        
        genconfdistsLT{igenmodel} = struct();
        genconfdistsTT{igenmodel} = struct();
        
        %%% GET CONFIDENCE COEFFICIENTS FROM DATA, MAKE GENERATIVE
        %%% DISTRIBUTIONS
        
        thiscoeffsLT = confcoeffLT{idcmodels(igenmodel)};
        thisRMSE = confRMSE{idcmodels(igenmodel)};
        
        %%% remove outliers
        idc = all(abs(zscore(thiscoeffsLT))<=1.96,2);
        thiscoeffsLT = thiscoeffsLT(idc,:);
        thisRMSE = thisRMSE(idc,:);
        %%%
        
        if doPostGen
            thiscoeffsTT = confcoeffTT{idcmodels(igenmodel)};
            thisRMSEPost = confRMSEPost{idcmodels(igenmodel)};
            
            %remove outliers
            idc = all(abs(zscore(thiscoeffsTT))<=1.96,2);
            thiscoeffsTT = thiscoeffsTT(idc,:);
            thisRMSEPost = thisRMSEPost(idc,:);
            %
        else
            thiscoeffsTT = zeros(size({idcmodels(igenmodel)}));
            thisRMSEPost = zeros(size(thisRMSE));
        end
        
        genconfdistsLT{igenmodel}.coeffs.intercept = @()random('Normal',mean(thiscoeffsLT(:,1)),std(thiscoeffsLT(:,1)));
        if doPostGen
            genconfdistsTT{igenmodel}.coeffs.intercept = @()random('Normal',mean(thiscoeffsTT(:,1)),std(thiscoeffsTT(:,1)));
        end
        %%% check which parameters to use based on coefficient names
        for iparam = 2:numel(coeffNames{igenmodel})
            genconfdistsLT{igenmodel}.coeffs.(coeffNames{igenmodel}{iparam}) =  @()random('Normal',mean(thiscoeffsLT(:,iparam)),std(thiscoeffsLT(:,iparam)));
            if doPostGen
                genconfdistsTT{igenmodel}.coeffs.(coeffNames{igenmodel}{iparam}) =  @()random('Normal',mean(thiscoeffsTT(:,iparam)),std(thiscoeffsTT(:,iparam)));
            end
        end
        
        pNL = lognfit(thisRMSE);
        pNT = lognfit(thisRMSEPost);
        
        genconfdistsLT{igenmodel}.confNoise = @()random('Lognormal',pNL(1),pNL(2));
        if doPostGen
            genconfdistsTT{igenmodel}.confNoise = @()random('Lognormal',pNT(1),pNT(2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        genconfnames = fieldnames(genconfdistsLT{igenmodel}.coeffs);
        
        for isub = 1:nsub
            paramstructgenconf =struct();
            paramstructgenconfpost =struct();
            for iparam = 1:numel(genconfnames)
                thisParam = genconfnames{iparam};
                paramstructgenconf.(thisParam) = genconfdistsLT{igenmodel}.coeffs.(thisParam)();
                genparamsconf{isim,isub,igenmodel} = [genparamsconf{isim,isub,igenmodel},paramstructgenconf.(thisParam)];
                if doPostGen
                    paramstructgenconfpost.(thisParam) = genconfdistsTT{igenmodel}.coeffs.(thisParam)();
                    genparamsconfpost{isim,isub,igenmodel} = [genparamsconfpost{isim,isub,igenmodel},paramstructgenconfpost.(thisParam)];
                end
            end
            %%% bookkeep noise
            gennoiseconf{isim,isub,igenmodel} = genconfdistsLT{igenmodel}.confNoise();
            gennoiseconfpost{isim,isub,igenmodel} = genconfdistsTT{igenmodel}.confNoise();
        end
    end
end

%% loop over simulations
%% generate RL time series for all subs x models

for isim = 1:nsims
    pars_sim = reshape(genparams(isim,:,whichlearnmodel),[nsub,1]);
    s_sim = reshape(s_all(isim,:,whichlearnmodel,:,:),[nsub,1,nsess,ntrials*4]);
    a_sim = reshape(choice_all(isim,:,whichlearnmodel,:,:),[nsub,1,nsess,ntrials*4]);
    r_sim = reshape(r_all(isim,:,whichlearnmodel,:,:),[nsub,1,nsess,ntrials*4]);
    c_sim = reshape(c_all(isim,:,whichlearnmodel,:,:),[nsub,1,nsess,ntrials*4]);
    ss_sim = reshape(ss_all(isim,:,whichlearnmodel,:,:),[nsub,1,2,ntrialspost]);
    aa_sim = reshape(choicePost_all(isim,:,whichlearnmodel,:),[nsub,1,ntrialspost]);
    
    %%% get RL vars
    [outRLSim,outRLSimSorted] = timeseriesRLPop(pars_sim,modelsinfo(whichlearnmodel),s_sim,a_sim,r_sim,c_sim,ss_sim,aa_sim);
    
    
    
    %%  simulate confidence for each model
    
    genparamsL = squeeze(genparamsconf(isim,:,:));
    genparamsT = squeeze(genparamsconfpost(isim,:,:));
    noiseL = squeeze(gennoiseconf(isim,:,:));
    noiseT = squeeze(gennoiseconfpost(isim,:,:));
    
    %%%transform input for simulateCOnfPop function
    %%% needs to be cell array of n models, each cell of size subj x pars
    
    for igenmodel = 1:nmodels
        genparsL2{igenmodel} = cell2mat(genparamsL(:,igenmodel));
        genparsT2{igenmodel} = cell2mat(genparamsT(:,igenmodel));
        noiseL2{igenmodel} = cell2mat(noiseL(:,igenmodel));
        noiseT2{igenmodel} = cell2mat(noiseL(:,igenmodel));
    end
    
    [confL,confT,confLsorted] = simulateConfPop(genparsL2,genparsT2,noiseL2,noiseT2,confmodelsinfo,confmodelsinfo,outRLSim);
    
    
    %%% make transfer matrices at least for debugging %%%
    ss_sub = permute(squeeze(outRLSim.ss),[1,3,2]);
    aa_sub = squeeze(outRLSim.choicePost);
    
    [prefMat_sim,accMat_sim,confChoiceMat_sim,confChoicePairMat_sim,pref_sim,accPost_sim,confPost_sim]= ...
        makePostMats(ss_sub,aa_sub,squeeze(confT(:,igenmodel,:)));
    
    %%%%store confidence %%%%
    conf_all(isim,:,:,:,:) = confL;
    confPost_all(isim,:,:,:) = confT;
    
    
    %% fit models
    [regL,regT] = fitConfModels(confmodelsinfo,outRLSim,confL,confT)
    
    %% store results of fitting on this generative model and simulation
    regress.BIC(isim,:,:,:) = regL.BIC;
    regress.BICPost(isim,:,:,:) = regT.BIC;
    regress.coeffs(isim,:,:,:) = regL.coef;
    regress.coeffsPost(isim,:,:,:) = regT.coef;
    regress.RMSE(isim,:,:,:) = regL.RMSE;
    regress.RMSEPost(isim,:,:,:) = regT.RMSE;
    
    if ~exist(['Results',filesep,outfilename,'.mat'],'file')
        save(['Results',filesep,outfilename,'.mat'], ...
            'genparamsconf','genparamsconfpost','regress','conf_all','confPost_all',...
            'modelsinfo','isim','-v7.3')
    else
        save(['Results',filesep,outfilename,'.mat'], ...
            'genparamsconf','genparamsconfpost','regress','conf_all','confPost_all',...
            'modelsinfo','isim','-append','-v7.3')
    end
    
end


function [regL,regT] = fitConfModels(confmodelsinfo,outRLSim,confL,confT)

nsub = size(confL,1);
nmodels = size(confmodelsinfo);

for igenmodel = 1:numel(confmodelsinfo)
    disp(['Gen Model ',num2str(igenmodel)])
    for isub = 1:size(confL,1)
        %%% LT vars and table %%%
        dQAbs_sub = abs(squeeze(outRLSim.dQ(isub,1,:,:,:)));
        Q_c_sub = squeeze(outRLSim.Qc(isub,1,:,:,:));
        Q_uc_sub = squeeze(outRLSim.Qu(isub,1,:,:,:));
        dQ_sub = Q_c_sub - Q_uc_sub;
        
        QdiffMaxMin_sub = squeeze(max(outRLSim.Q(isub,1,:,:,:),[],4)- min(outRLSim.Q(isub,1,:,:,:,:),[],4));
        
        V_sub = squeeze(outRLSim.V(isub,1,:,:,:));
        QdiffGB_sub =squeeze(outRLSim.Q(isub,1,:,2,:)-outRLSim.Q(isub,1,:,1,:)); %good minus bad
        pc_sub = squeeze(outRLSim.pc(isub,:,:,:));
        
        conf_mat_sub = squeeze(confL(isub,igenmodel,:,:));
        
        conf_prev_mat = [0,0,0;conf_mat_sub(:,1:end-1)']';
        
        conf_mat_sub(conf_mat_sub==1) = 0.9999;
        conf_prev_mat(conf_prev_mat==1) = 0.9999;
        conf_prev_mat(conf_prev_mat==0) = 1-0.9999;
        
        varnames = {'conf','dQ','dQabs','Qc','Qu','sigmaQ','V','confPrev','dQGoodBad','dQMaxMin','pc'};
        tblL = table(conf_mat_sub(:),dQ_sub(:),abs(dQ_sub(:)),Q_c_sub(:),Q_uc_sub(:),...
            [Q_c_sub(:)+Q_uc_sub(:)],V_sub(:),conf_prev_mat(:),QdiffGB_sub(:),QdiffMaxMin_sub(:),pc_sub(:),...
            'VariableNames', varnames);
        %         tbl = table(zscore(conf_mat_sub(:)),zscore(dQ_sub(:)),zscore(abs(dQ_sub(:))),zscore(V_sub(:)),zscore(Q_c_sub(:)),zscore(Q_uc_sub(:)),zscore([Q_c_sub(:)+Q_uc_sub(:)]),zscore(conf_prev_mat(:)),'VariableNames', varnames);
        
        %%% TT vars and table %%%%%
        dQ_post_sub = squeeze(outRLSim.dQ_post(isub,1,:));
        Q_c_post_sub = squeeze(outRLSim.Qc_post(isub,1,:));
        Q_uc_post_sub = squeeze(outRLSim.Qu_post(isub,1,:));
        V_post_sub = squeeze(outRLSim.V_post(isub,1,:));
        dQGoodBad_post_sub = squeeze(outRLSim.Q_post(isub,1,2,:)-outRLSim.Q_post(isub,1,1,:));
        pc_post_sub = squeeze(outRLSim.pc_post(isub,1,:));
        
        sigmaQ_post = Q_c_post_sub + Q_uc_post_sub;
        
        confPost_sub = squeeze(confT(isub,igenmodel,:));
        confPost_sub(confPost_sub==1) = 0.9999;
        confPostPrev = [0,confPost_sub(1:end-1)']';
        confPostPrev(confPostPrev==0) = 1-0.9999;
        
        varnames = {'conf','dQ','dQabs','Qc','Qu', 'sigmaQ','V','confPrev','dQGoodBad','pc'};
        tblT = table([confPost_sub(:)],[dQ_post_sub(:)],[abs(dQ_post_sub(:))],[Q_c_post_sub(:)],...
            [Q_uc_post_sub(:)],[Q_c_post_sub(:)]+[Q_uc_post_sub(:)],V_post_sub(:),confPostPrev(:),...
            dQGoodBad_post_sub(:),pc_post_sub(:),'VariableNames', varnames);
        
        
        %%% fit all conf models %%%
        parfor irecmodel =1:numel(confmodelsinfo)
            %             disp(['recmodel ', num2str(irecmodel)]);
            formulaCONF = 'conf ~ 1';
            for ipar=2:numel(confmodelsinfo{irecmodel}.varnames)
                formulaCONF = [formulaCONF,'+',confmodelsinfo{irecmodel}.varnames{ipar}];
            end
            
            glm = fitglm(tblL,formulaCONF,'link','logit');
            AICL(isub,igenmodel,irecmodel) =  glm.ModelCriterion.AIC;
            BICL(isub,igenmodel,irecmodel) =  glm.ModelCriterion.BIC;
            coefL{isub,igenmodel,irecmodel} = glm.Coefficients.Estimate;
            RMSEL(isub,igenmodel,irecmodel) = nanstd(glm.Residuals.LinearPredictor);
            
            %%%% Transfer
            glm = fitglm(tblT,formulaCONF,'link','logit');
            AICT(isub,igenmodel,irecmodel) =  glm.ModelCriterion.AIC;
            BICT(isub,igenmodel,irecmodel) =  glm.ModelCriterion.BIC;
            coefT{isub,igenmodel,irecmodel} = glm.Coefficients.Estimate;
            RMSET(isub,igenmodel,irecmodel) = nanstd(glm.Residuals.LinearPredictor);
        end
    end
end

regL.BIC=BICL;regL.AIC =AICL; regL.coef = coefL;regL.RMSE = RMSEL;
regT.BIC=BICT;regT.AIC =AICT; regT.coef = coefT;regT.RMSE = RMSET;

end

