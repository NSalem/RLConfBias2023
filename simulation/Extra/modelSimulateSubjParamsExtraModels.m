%%% simulate choice and confidence based on subject parameters %%

addpath('helperfuncs')
addpath('ModelingFuncs')

% load learning parameters
load('Results\RLvars_ExtraContext_all.mat', 'params_allexp','modelsinfo')



rlnames = {'ROther2','ROtherNorm','Perseveration'};
for learnmodel = 1:3
    
    % learnmodel = 3;
    
    reg = load(['Results\reg_conflogit',rlnames{learnmodel},'_learning_dqabs.mat']);
    regPost = load(['Results\reg_conflogit',rlnames{learnmodel},'_posttest_dqabs.mat']);
    idcLT = 7:10;%find(reg.whichLearnModel ==learnmodel & reg.whichConfModel>1 & reg.isConfPrev); %select conf model
    idcTT = 1:4;%find(regPost.whichLearnModel ==learnmodel & regPost.whichConfModel>1 & ~regPost.isConfPrev); %select conf model
    
    %%% use subject parameters
    learnparams = cell2mat({params_allexp{:,learnmodel}}'); %select learn model
    
    nsims = 20;
    
    %%%%% subject parameters
    confparamsLT = cell([numel(idcLT),1]);
    confNoiseLT = cell([numel(idcLT),1]);
    nsub = size(learnparams,1);
    
    for iconfmodelLT = 1:numel(idcLT)
        confparamsLT{iconfmodelLT} = reg.confcoeff{idcLT(iconfmodelLT)};
        confNoiseLT{iconfmodelLT} =reg.confRMSE{idcLT(iconfmodelLT)};
    end
    
    confparamsTT = cell([numel(idcTT),1]);
    confNoiseTT = cell([numel(idcTT),1]);
    
    for iconfmodelTT = 1:numel(idcTT)
        confparamsTT{iconfmodelTT} =  regPost.confcoeff{idcTT(iconfmodelTT)};
        confNoiseTT{iconfmodelTT} = regPost.confRMSE{idcTT(iconfmodelTT)};
    end
    
    %% sample rl params from fitted distribution
    nsub = size(learnparams,1);
    nsub = 1000;
    
    learnparamssub = learnparams;
    learnparams = nan(nsub,size(learnparams,2));
    for iparam = 1:size(learnparams,2)
        if iparam ==1
            p = gamfit(learnparamssub(:,iparam));
            learnparams(:,iparam) = gamrnd(p(1),p(2),size(learnparams,1),1);
        elseif iparam == 6
            p(1),p(2) = normfit(learnparamssub(:,iparam));
            learnparams(:,iparam) = betarnd(p(1),p(2),size(learnparams,1),1);
        else
            p = betafit(learnparamssub(:,iparam));
            learnparams(:,iparam) = betarnd(p(1),p(2),size(learnparams,1),1);
        end
    end
    
    %%%%%%
    learncell = cell([size(learnparams,1),1]);
    for i = 1:numel(learncell)
        learncell{i,1} = learnparams(i,:);
    end
    
    %% random conf params with subj gaussian distributions
    confparamsLT = cell([numel(idcLT),1]);
    confNoiseLT = cell([numel(idcLT),1]);
    
    for iconfmodelLT = 1:numel(idcLT)
        
        coeffsThisModel = reg.confcoeff{idcLT(iconfmodelLT)};
        idc = all(abs(zscore(coeffsThisModel))<=1.96,2);
        pm = mean(coeffsThisModel(idc,:));
        ps = std(coeffsThisModel(idc,:));
        
        confparamsLT{iconfmodelLT} =  randn(nsub,size(reg.confcoeff{idcLT(iconfmodelLT)},2)).*ps+pm;
        
        p = lognfit(reg.confRMSE{idcLT(iconfmodelLT)}(idc,:));
        confNoiseLT{iconfmodelLT} =lognrnd(p(1),p(2),nsub,1);
    end
    
    confparamsTT = cell([numel(idcTT),1]);
    confNoiseTT = cell([numel(idcTT),1]);
    
    for iconfmodelTT = 1:numel(idcTT)
        
        %%% remove outliers
        coeffsThisModel = regPost.confcoeff{idcTT(iconfmodelTT)};
        idc = all(abs(zscore(coeffsThisModel))<=1.96,2);
        pm = mean(coeffsThisModel(idc,:));
        ps = std(coeffsThisModel(idc,:));
        
        confparamsTT{iconfmodelTT} =  randn(nsub,size(regPost.confcoeff{idcTT(iconfmodelTT)},2)).*ps+pm;
        
        p = lognfit(regPost.confRMSE{idcTT(iconfmodelTT)}(idc,:));
        confNoiseTT{iconfmodelTT} =lognrnd(p(1),p(2),nsub,1);
    end
    
    %%%
    outcomes = [0.1,1];
    ntrials = 24;
    nsess = 3;
    paramstruct = modelsinfo{learnmodel};
    
    optP = [0.75,0.25,0.75,0.25,0.25,0.75,0.25,0.75];
    optU = [optP*1+(1-optP)*0.1]; %option utility (values in exp 1 are different but same order)
    optU(5:8) = -optU(5:8);
    % OUT = Generate_Outcomes(ntrials,outcomes);
    
    %%% do simulations
    %%for isim = 1:50
    
    %generate parameters
    % genpars = getlearnparams(modelsinfo,10,
    
    
    confmodelsinfo = cell(1,1);
    confmodelsinfo{1} = struct()
    
    cfgL  = cell(1,1);
    cfgT = cell(1,1);
    
    confmodelsinfo{1} = struct('varnames',{{'intercept','dQabs','confPrev'}},'ilearnmodel',1);
    confmodelsinfo{2} = struct('varnames',{{'intercept','dQabs','Qc','confPrev'}},'ilearnmodel',1);
    confmodelsinfo{3} = struct('varnames',{{'intercept','dQabs','sigmaQ','confPrev'}},'ilearnmodel',1);
    confmodelsinfo{4} = struct('varnames',{{'intercept','dQabs','V','confPrev'}},'ilearnmodel',1);
    
    
    confmodelsinfoT{1} = struct('varnames',{{'intercept','dQabs'}},'ilearnmodel',1);
    confmodelsinfoT{2} = struct('varnames',{{'intercept','dQabs','Qc'}},'ilearnmodel',1);
    confmodelsinfoT{3} = struct('varnames',{{'intercept','dQabs','sigmaQ'}},'ilearnmodel',1);
    confmodelsinfoT{4} = struct('varnames',{{'intercept','dQabs','V'}},'ilearnmodel',1);
    
    
    %%% loop for n simulations
    outConfL = nan([nsims,nsub,numel(idcLT),3,96]);
    outConfT = nan([nsims,nsub,numel(idcTT),112]);
    outConfLSorted = nan([nsims,nsub,numel(idcLT),3,4,24]);
    prefMat = nan([nsims,nsub,8,8]);
    accMat = nan([nsims,nsub,8,8]);
    confChoiceMat = nan([nsims,numel(idcTT),nsub,8,8]);
    confChoicePairMat = nan([nsims,numel(idcTT),nsub,8,8]);
    %
    
    for isim = 1:nsims
        
        %simulate RL
        [outRL,outSorted]= simulateRLPop(learncell,modelsinfo(learnmodel),3, 24, outcomes);
        
        %simulate confidence
        
        [outConfL_sim,outConfT_sim,outConfLSorted_sim] = simulateConfPop(confparamsLT,confparamsTT,confNoiseLT,confNoiseTT,...
            confmodelsinfo,confmodelsinfoT,outRL);
        
        %%%  make post mats for preference, accuracy and for each conf model %%%
        ss_sub = permute(squeeze(outRL.ss),[1,3,2]);
        aa_sub = squeeze(outRL.choicePost);
        for imodel = 1:numel(idcTT)
            conf_sub = squeeze(outConfT_sim);
            %         [prefMat(isim,:,:,:),accMat(isim,:,:,:),confChoiceMat(isim,imodel,:,:,:),confChoicePairMat(isim,imodel,:,:,:)]= makePostMats(ss_sub,aa_sub,squeeze(outConfT_sim(:,imodel,:)));
            
            [prefMat_sim,accMat_sim,confChoiceMat_sim,confChoicePairMat_sim,pref_sim,accPost_sim,confPost_sim]= ...
                makePostMats(ss_sub,aa_sub,squeeze(outConfT_sim(:,imodel,:)));
            
            prefMat(isim,:,:,:) = prefMat_sim;
            accMat(isim,:,:,:) = accMat_sim;
            confChoiceMat(isim,imodel,:,:,:) = confChoiceMat_sim;
            confChoicePairMat(isim,imodel,:,:,:) = confChoicePairMat_sim;
            pref(isim,:,:) = pref_sim;
            accPost(isim,:,:) = accPost_sim;
            confPost(isim,imodel,:,:) = confPost_sim;
            
        end
        %%%
        
        outConfL(isim,:,:,:,:) = outConfL_sim;
        outConfT(isim,:,:,:) = outConfT_sim;
        outConfLSorted(isim,:,:,:,:,:) = outConfLSorted_sim;
    end
    save(['Results\SimulationsSubjParams',rlnames{learnmodel},'.mat'],'outRL','outSorted',...
        'outConfL','outConfT','outConfLSorted','modelsinfo','learnparams','confparamsLT','confNoiseLT','confparamsTT','confNoiseTT',...
        'prefMat','accMat','confChoiceMat','confChoicePairMat','pref','accPost','confPost');
end