function regressConfLT(rlvarspath,learnmodels,confmodels,outpath)
% fits regression models on confidence for the Learning Task and save to 
% specified path
% inputs:
%- rlvarspath: RLVars file path
%- learnmodels: list of learning models to use
%- confmodels: list of confidence models,as the right hand side of regression formula 
%- outpath: where to save results

load(rlvarspath);
nsub = size(dQ,2);
% nmodels = size(
nsess = size(dQ,3);
ntrials = size(dQ,5);
confmodel = nan([nsub,10,nsess,4,ntrials]);
confcoeff = cell(10,1);
conft = cell(10,1);
confRMSE = cell(10,1);

BICCONF = []; %initialize var

for isub = 1:size(correct,1)

    disp(['Participant ',num2str(isub)]);
    imodelcomb = 0; %

    for ilearnmodel = learnmodels
        %% Get hidden variables for this model for this subject
%         dQ_sub = squeeze(dQ(ilearnmodel,isub,:,:,:));
        dQAbs_sub = abs(squeeze(dQ(ilearnmodel,isub,:,:,:)));
        Q_c_sub = squeeze(Q_c(ilearnmodel,isub,:,:,:));
        Q_uc_sub = squeeze(Q_uc(ilearnmodel,isub,:,:,:));
        dQ_sub = Q_c_sub - Q_uc_sub;      
        QdiffMaxMin_sub = squeeze(max(Q(ilearnmodel,isub,:,:,:,:),[],5)- min(Q(ilearnmodel,isub,:,:,:,:),[],5));      
        V_sub = squeeze(V(ilearnmodel,isub,:,:,:));
        QdiffGB_sub =squeeze(Q(ilearnmodel,isub,:,:,2,:)-Q(ilearnmodel,isub,:,:,1,:)); %good minus bad
        pc_sub = squeeze(pc(ilearnmodel,isub,:,:,:));
             
        %% Regressions  
        conf_mat_sub = squeeze(conf(isub,:,:,:));
        conf_prev_mat = squeeze(confprev(isub,:,:,:));
        
        conf_mat_sub2 = (conf_mat_sub-0.5)*2;
        conf_prev_mat2 = (conf_prev_mat-0.5)*2;

        conf_mat_sub(conf_mat_sub==1) = 0.9999;
        conf_mat_sub(conf_mat_sub==0) = 1-0.9999;
        conf_prev_mat(conf_prev_mat==1) = 0.9999;
        conf_prev_mat(conf_prev_mat==0) = 1-0.9999;
        
        conf_mat_sub2(conf_mat_sub2==1) = 0.9999;
        conf_mat_sub2(conf_mat_sub2==0) = 1-0.9999;
        conf_prev_mat2(conf_prev_mat2==1) = 0.9999;
        conf_prev_mat2(conf_prev_mat2==0) = 1-0.9999;

        varnames = {'conf','conf2','dQ','dQabs','Qc','Qu','QcplusQu','V','confprev','dQGoodBad','dQMaxMin','pc'};
        tbl = table(conf_mat_sub(:),conf_mat_sub2(:),dQ_sub(:),abs(dQ_sub(:)),Q_c_sub(:),Q_uc_sub(:),...
        [Q_c_sub(:)+Q_uc_sub(:)],V_sub(:),conf_prev_mat(:),QdiffGB_sub(:),QdiffMaxMin_sub(:),pc_sub(:),...
        'VariableNames', varnames);

        for iconfmodel =1:numel(confmodels)
            imodelcomb =imodelcomb+1;
            
            formulaCONF = ['conf ~ 1 + ',confmodels{iconfmodel}];

            formulas{imodelcomb}  = formulaCONF;
            glm = fitglm(tbl,formulaCONF,'link','logit');     

            glmCONF{isub, imodelcomb} =glm;
            VIF{isub,imodelcomb} = diag(inv(corrcov(glm.CoefficientCovariance)));
            confmodel(isub,imodelcomb,:,:,:) = reshape(glm.predict(tbl),[nsess,4,ntrials]);
            BICCONF(isub,imodelcomb) = glm.ModelCriterion.BIC;
            LLCONF(isub,imodelcomb) = glm.LogLikelihood;
            npars(isub,imodelcomb) = glm.NumCoefficients;
            confcoeff{imodelcomb}(isub,:) = glm.Coefficients.Estimate;
            
            conft{imodelcomb}(isub,:) = glm.Coefficients.tStat;
            confRMSE{imodelcomb}(isub,:) = nanstd(glm.Residuals.LinearPredictor);
            
            %%% also linear 
            lm = fitlm(tbl,formulaCONF);
            lmCONF{isub, imodelcomb} =lm;
            confmodel_linear(isub,imodelcomb,:,:,:) = reshape(lm.predict(tbl),[nsess,4,ntrials]);
            BICCONF_linear(isub,imodelcomb) = lm.ModelCriterion.BIC;

            confcoeff_linear{imodelcomb}(isub,:) = lm.Coefficients.Estimate;
            conft_linear{imodelcomb}(isub,:) = lm.Coefficients.tStat;
            confRMSE_linear{imodelcomb}(isub,:) = lm.RMSE;
          
            %%% bookkeeping of model specifications
            modelsinfo = infostr.modelsinfo{infostr.subexp(isub)};
            isLearnRel(imodelcomb) = isfield(modelsinfo{ilearnmodel},'contextual')&&modelsinfo{ilearnmodel}.contextual;
            %                         isConfBiased(imodelcomb) = iregcontext>1;
            whichConfModel(imodelcomb) = iconfmodel;
            confBias{imodelcomb} = erase(formulaCONF,{'conf ~ 1','+ dQabs','+ dQGoodBad', '+ dQMaxMin','+ dQ','+ pc','+ confprev',' '}); 
            whichLearnModel(imodelcomb) = ilearnmodel;
            isConfPrev(imodelcomb) = contains(formulaCONF,'confprev');
        end
    end
end

clear V_sub dQ_sub Q_c_sub Q_uc_sub
confmat = conf;
fitinfo.date = date;
save(outpath,'glmCONF','BICCONF','LLCONF','confcoeff','conft','isLearnRel','confBias','whichConfModel','whichLearnModel','isConfPrev','confmat', 'confmodel','confRMSE','fitinfo','formulas','VIF','npars','BICCONF_linear','confcoeff_linear','conft_linear','confRMSE_linear')
clear formulas glm glmCONF confmodel BICCONF  LLCONF confcoeff conft conftconfRMSE islearnRel isConfBias whichConfModel whichLearnModel isConfPrev
end
