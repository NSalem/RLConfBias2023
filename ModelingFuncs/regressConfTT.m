function regressConfTT(rlvarspath,learnmodels,confmodels,outpath)
% fits regression models on confidence for the Transfer Task and save to 
% specified path
% inputs:
%- rlvarspath: RLVars file path
%- learnmodels: list of learning models to use
%- confmodels: list of confidence models,as the right hand side of regression formula 
%- outpath: where to save results
load(rlvarspath);

% nsub = size(dQ,2);

BICCONF = []; %initialize var

for isub = 1:size(correct,1)
    
    disp(['Participant ',num2str(isub)]);
    imodelcomb = 0;
    
    %% load actual confidence
    conf = confpost(isub,:);
    s1 = size(conf,1);
    s2 = size(conf,2);
    confprev = [0,conf(:,1:end-1)];
    
    conf(conf==1) = 0.9999;
    confprev(confprev==1) = 0.9999;
    confprev(confprev==0) = 1-0.9999;
    
    for ilearnmodel = learnmodels
        ss_sub =squeeze(ss(isub,:,:));
        aa_sub = squeeze(aa(isub,:));
        %% Calculate hidden variables
        dQ_post_sub = squeeze(dQ_post(ilearnmodel,isub,:));
        Q_c_post_sub = squeeze(Q_c_post(ilearnmodel,isub,:));
        Q_uc_post_sub = squeeze(Q_uc_post(ilearnmodel,isub,:));
        V_post_sub = squeeze(V_post(ilearnmodel,isub,:));
        dQGoodBad_sub = squeeze(QG_post(ilearnmodel,isub,:)-QB_post(ilearnmodel,isub,:));
        pc_post_sub = squeeze(pc_post(ilearnmodel,isub,:));
        %% Regressions
        varnames = {'conf','dQ','dQabs','Qc','Qu', 'QcplusQu','V','confprev','dQGoodBad','pc'};
        tbl = table([conf(:)],[dQ_post_sub(:)],[abs(dQ_post_sub(:))],[Q_c_post_sub(:)],...
            [Q_uc_post_sub(:)],[Q_c_post_sub(:)]+[Q_uc_post_sub(:)],V_post_sub(:),confprev(:),...
             dQGoodBad_sub,pc_post_sub,'VariableNames', varnames);
        
        for iconfmodel =1:numel(confmodels)
            imodelcomb =imodelcomb+1;
            
            formulaCONF = ['conf ~ 1 + ',confmodels{iconfmodel}];
            formulas{imodelcomb} = formulaCONF;
            thisReg = fitglm(tbl,formulaCONF,'link','logit');
            VIF{isub,imodelcomb} = diag(inv(corrcov(thisReg.CoefficientCovariance)));
            confBias{imodelcomb} = erase(formulaCONF,{'conf ~ 1','+ dQabs','+ dQGoodBad','+ dQ','+ pc','+ confprev',' '});
            confmodel(isub,imodelcomb,:) = thisReg.predict(tbl);
            glmCONF{isub, imodelcomb} = thisReg;

            BICCONF(isub,imodelcomb) = thisReg.ModelCriterion.BIC;
            LLCONF(isub,imodelcomb) = thisReg.LogLikelihood;
            npars(isub,imodelcomb) = thisReg.NumCoefficients;
            confcoeff{imodelcomb}(isub,:) = thisReg.Coefficients.Estimate;
            confRMSE{imodelcomb}(isub,:) =   nanstd(thisReg.Residuals.LinearPredictor);
            conft{imodelcomb}(isub,:) = thisReg.Coefficients.tStat;
            
            
            %%% also linear
            thisReg = fitlm(tbl,formulaCONF);
            confmodel_linear(isub,imodelcomb,:) = thisReg.predict(tbl);
            
            lmCONF{isub, imodelcomb} = thisReg;
            BICCONF_linear(isub,imodelcomb) = thisReg.ModelCriterion.BIC;
            confcoeff_linear{imodelcomb}(isub,:) = thisReg.Coefficients.Estimate;
            confRMSE_linear{imodelcomb}(isub,:) =   nanstd(thisReg.RMSE);
            conft_linear{imodelcomb}(isub,:) = thisReg.Coefficients.tStat;
   
            %%%
            
            %%% bookkeeping 
            modelsinfo = infostr.modelsinfo{infostr.subexp(isub)};
            isLearnRel(imodelcomb) = isfield(modelsinfo{ilearnmodel},'contextual')&&modelsinfo{ilearnmodel}.contextual;
%             isConfBiased(imodelcomb) = iregcontext>1;
            whichConfModel(imodelcomb) = iconfmodel;
            whichLearnModel(imodelcomb) = ilearnmodel;
            isConfPrev(imodelcomb) = contains(formulaCONF,'confprev');
        end
        
        Qc_model(isub,ilearnmodel,:) = Q_c_post_sub;
        dQ_model(isub,ilearnmodel,:) = dQ_post_sub;
        
    end
    
    clear V dQ Q_c Q_uc
    
    %%  make matrix of confidence for each pair of symbols 
    for k_symb = 1:8
        for k_symb2 = 1:8
            if k_symb == k_symb2
                confmat(isub,k_symb2,k_symb) = NaN;
                pref(isub,k_symb2,k_symb)= NaN;
                Ppref(isub,k_symb2,k_symb) = NaN;
                Ppref2(isub,k_symb2,k_symb) = NaN;
                dQmat(isub,k_symb2,k_symb) = NaN;
                Vmodel(isub,k_symb2,k_symb) = NaN;
                Qcmat(isub,k_symb2,k_symb) = NaN;
            else
                S1 = (ss_sub(:,1) == k_symb) & (ss_sub(:,2) == k_symb2);
                S2 = (ss_sub(:,2) == k_symb) & (ss_sub(:,1) == k_symb2);
                
                confmat(isub,k_symb2,k_symb) = (nansum(conf(S1)) + nansum(conf(S2)))./(4); %actual confidence
                confmatChoice(isub,k_symb2,k_symb) = nanmean(conf((S1& aa_sub' ==1 )|(S2&aa_sub' ==2)));
%                 confChoiceMat(isub,k_symb2,k_symb) = (nanmean(conf(isub,(S1& aa_sub' ==1 )|(S2&aa_sub' ==2))));            

                for imodel = 1:imodelcomb %confidence for each model
                    confmodelmat(imodel,isub,k_symb2,k_symb) = (nansum(confmodel(isub,imodel,S1)) + nansum(confmodel(isub,imodel,S2)))./(4);
                    
                    confmodelmatChoice(imodel,isub,k_symb2,k_symb) = (nanmean(confmodel(isub,imodel,(S1& aa_sub' ==1 )|(S2&aa_sub' ==2))));
                    confmodelmatChoice(imodel,isub,k_symb2,k_symb) = (nanmean(confmodel(isub,imodel,(S1& aa_sub' ==1 )|(S2&aa_sub' ==2))));
                    %                             dumDQ(imodel,kk,k_symb2,k_symb) =  (nanmean(abs(dQ_post_sub((S1& aa' ==1 )|(S2&aa' ==2)))));
                    %                             if nanmean(abs(dQ_post_sub((S1& aa' ==1 )|(S2&aa' ==2))))~=nanmean(abs(dQ_post_sub((S1)|(S2))))
                    %                             keyboard()
                    %                             end
                end
                
                dQmat(isub,k_symb2,k_symb) = nansum((dQ_model(isub,ilearnmodel,S2))) ./(2);
                Qcmat(isub,k_symb2,k_symb) = (nansum(Qc_model(isub,ilearnmodel,S2)))./(2);
                
            end
        end
    end
    
    %% confidence when choosing each symbol
    for k_symb = 1:8;
        S1 = ss_sub(:,1) == k_symb;
        S2 = ss_sub(:,2)== k_symb;
        %                 pref(kk,k_symb) = 100*(mean(aa(S1)==1) + mean(aa(S2)==2))./2;
        confpost(isub,k_symb) =squeeze(nanmean(conf((S1 & aa_sub' ==1) | (S2 & aa_sub' ==2)))); %actual confidence
        confmodelpost(isub,:,k_symb) = squeeze(nanmean(confmodel(isub,:,(S1 & aa_sub' ==1) | (S2 & aa_sub' ==2)),3));
        confmodelpostPresent(isub,:,k_symb) = squeeze(nanmean(confmodel(isub,:,(S1|S2)),3));
    end
    
end
clear parametersLPP
save(outpath,'BICCONF','LLCONF','confcoeff','conft','confRMSE','isLearnRel','whichConfModel','whichLearnModel','isConfPrev','confpost','confmat', 'confmatChoice','confmodelpost','confmodelpostPresent','confmodel','confmodelmatChoice','confmodelmat','formulas','confBias','VIF','BICCONF_linear','confcoeff_linear','conft_linear','confRMSE_linear','npars')
clear formulas confmodel confmodelDiscrete confmodelmat confmodelmatChoice confmat pref dQmat Vmodel Qcmat Qc_model confpost confmodelpost confmodelpost BICCONF LLCONF confcoeff confRMSE conft isLearnRel isConfPrev whichConfModel whichLearnModel isConfBiased
end
