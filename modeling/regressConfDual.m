
addpath('ModelingFuncs\');
addpath('helperfuncs');
resultsdir = ['Results',filesep];
datadir = ['data',filesep];
outfilename = '';

ilearnmodel = 11; 
whichexp = 1:5; % select experiments


exps = load([resultsdir,'data_all_CONF']);
load([resultsdir,'RLvars']); %load hidden variables from learning models
exps_conf = exps.data_all;


BICCONF = []; %initialize var

kk = 0;
for iexp = whichexp
    thisExp = MakeMatrix_sess_cond_trl(exps_conf(iexp),1); %make matrices of subj x sess x cond x trial
    modelsinfo = infostr.modelsinfo{iexp};    
    for isub = 1:size(exps_conf(iexp).corr,1) 
        kk = kk+1; % subject number accross all experiments
        imodelcomb = 0; %
        
        conf_mat = squeeze(thisExp.conf_mat_sess(isub,:,:,:));
        rt_mat   = squeeze(thisExp.RT_mat(isub,:,:,:));
        conf_prev_mat = squeeze(thisExp.prevconf_mat_sess(isub,:,:,:));
        rtPost = exps_conf(iexp).rt_post(isub,:);
        confPost = exps_conf(iexp).conf_post(isub,:);
        confprevPost = [0,confPost(:,1:end-1)];
        
        conf_mat(conf_mat==1) = 0.9999;
        conf_mat(conf_mat==0) = 1-0.9999;
        conf_prev_mat(conf_prev_mat==1) = 0.9999;
        conf_prev_mat(conf_prev_mat==0) = 1-0.9999;
       
                
        confPost(confPost==1) = 0.9999;
        confPost(confPost==0) = 1-0.9999;
        confprevPost(confprevPost==1) = 0.9999;
        confprevPost(confprevPost==0) = 1-0.9999;
       
        
        %% Get hidden variables for this model for this subject 
        dQ_sub = squeeze(dQ(ilearnmodel,kk,:,:,:));
        Q_c_sub = squeeze(Q_c(ilearnmodel,kk,:,:,:));
        Q_uc_sub = squeeze(Q_uc(ilearnmodel,kk,:,:,:));
        V_sub = squeeze(V(ilearnmodel,kk,:,:,:));

        ss = squeeze(exps_conf(iexp).ss(isub,:,:));
        aa = squeeze(exps_conf(iexp).aa(isub,:));
        %% Calculate hidden variables for transfer task
        dQ_post_sub = squeeze(dQ_post(ilearnmodel,kk,:));
        V_post_sub = squeeze(V_post(ilearnmodel,kk,:));
        Q_c_post_sub = squeeze(Q_c_post(ilearnmodel,kk,:));
        Q_uc_post_sub = squeeze(Q_uc_post(ilearnmodel,kk,:));
        
        nL = numel(dQ_sub(:)); %ntrials LT;
        nT = numel(dQ_post_sub(:));%ntrials TT;
        varnames = {'conf','rt','intercept','dQ','Qc','sigmaQ','V','intercept_L','dQ_L','Qc_L','sigmaQ_L','V_L','intercept_T','dQ_T','Qc_T','sigmaQ_T','V_T','confprev'};
        tbl = table([conf_mat(:);confPost'],...%conf
            [rt_mat(:);rtPost'],... %RT
            [ones(nL,1);ones(numel(confPost),1)],... %intercept
            abs([dQ_sub(:);dQ_post_sub(:)]),... %dQ
            [Q_c_sub(:);Q_c_post_sub(:)],... %Qc
            [Q_c_sub(:)+Q_uc_sub(:);Q_c_post_sub(:)+Q_uc_post_sub(:)],...%sigmaQ
            [V_sub(:);V_post_sub(:)],... %V
            [ones(nL,1);zeros(numel(confPost),1)],... %intercept_L            
            abs([dQ_sub(:);zeros(nT,1)]),... dQ_L
            [Q_c_sub(:);zeros(nT,1)],... %Qc_L
            [Q_uc_sub(:)+Q_c_sub(:);zeros(nT,1)],...%sigmaQ_L
            [V_sub(:);zeros(nT,1)],... %V_L
            [zeros(nL,1);ones(numel(confPost),1)],... %intercept_T
            abs([zeros(nL,1);dQ_post_sub(:)]),... %dQ_T
            [zeros(nL,1);Q_c_post_sub(:)],... %Qc_T
            [zeros(nL,1);Q_uc_post_sub(:)+Q_c_post_sub(:)],... %sigmaQ_T
            [zeros(nL,1);V_post_sub(:)],... %V_T
            [conf_prev_mat(:);zeros(nT,1)],... %confprev only for LT, given BMC results
            'VariableNames', varnames);
            
        
        imodelcomb = 0; %combined model index, total 32
        for idiffB0 = 1:2 % 2: different intercept for transfer task 
            if idiffB0 == 2
                offsetReg = 'intercept_L + intercept_T';
            else 
                offsetReg = 'intercept';
            end
            for idiffBDQ = 1:2 %2: different Beta_dQ for transfer task 
                if idiffBDQ == 2
                    diffReg = 'dQ_L + dQ_T';
                else
                    diffReg = 'dQ';
                end
                for iregcontextLT = 1:3 %context bias for LT. 1:Qc,2:sigmaQ,3:V
                    for iregcontextTT = 1:3 %context bias for TT. 1:Qc, 2:sigmaQ                  
                        for idiffpar = 1:2 %2: using different parameters (if iregcontextLT = iregcontextTT)
                            switch iregcontextLT % determine context regressor for LT
                                case 1
                                    contextRegLT = 'Qc';
                                case 2
                                    contextRegLT = 'sigmaQ';
                                case 3
                                    contextRegLT = 'V';
                            end
                            
                            switch iregcontextTT %determine context regressor for TT
                                case 1
                                    contextRegTT = 'Qc';
                                case 2
                                    contextRegTT = 'sigmaQ';
                                case 3
                                    contextRegTT = 'V';
                            end
                            % decide general context regressor (whole task)
                            % vs per-task context regressors
                            if (idiffpar==2 && iregcontextLT == iregcontextTT) || (idiffpar == 1 && iregcontextLT ~= iregcontextTT)                   
                                if ~isempty(contextRegLT) %&& ~strcmp(contextRegLT,'V')
                                    contextRegLT = [contextRegLT,'_L'];
                                end
                                if ~isempty(contextRegTT) %&& ~strcmp(contextRegTT,'V')
                                    contextRegTT = [contextRegTT,'_T'];
                                end
                                contextReg = [contextRegLT,'+', contextRegTT];
                            elseif idiffpar == 1 && iregcontextLT == iregcontextTT
                                contextReg = contextRegLT;
                            end
                            formulaCONF = char(['conf~',offsetReg,'+',diffReg,'+',contextReg,'+ confprev']);
                            formulaRT = char(['conf~',offsetReg,'+',diffReg,'+',contextReg]);
 
                            if (idiffpar==2 && iregcontextLT == iregcontextTT) || idiffpar ==1 % idffpart == 2 only applies when using the same context model for LT and TT) 
                                imodelcomb = imodelcomb+1;        
                                
                                thisRegGlm = fitglm(tbl,formulaCONF,'Intercept',false,'link','logit');
                                thisRegLin= fitlm(tbl,formulaCONF,'Intercept',false);
                                
                                confmodelLT(kk,imodelcomb,:,:,:) = reshape(thisRegGlm.predict(tbl(1:288,:)),[3,4,24]);
                                confmodelTT(kk,imodelcomb,:) = thisRegGlm.predict(tbl(289:end,:));
                                glmCONF{iexp}{isub, imodelcomb} =thisRegGlm;
                                lmCONF{iexp}{isub, imodelcomb} =thisRegLin;
                                AICCONF(kk,imodelcomb) = thisRegGlm.ModelCriterion.AIC;
                                BICCONF(kk,imodelcomb) = thisRegGlm.ModelCriterion.BIC;
                                confcoeff{imodelcomb}(kk,:) = thisRegGlm.Coefficients.Estimate;
                                conft{imodelcomb}(kk,:) = thisRegGlm.Coefficients.tStat;
                                confRMSE{imodelcomb}(kk,:) =  nanstd(thisRegGlm.Residuals.LinearPredictor);

%                                 thisRegGlm =  fitglm(tbl,formulaCONF,'Intercept',false,'link','logit');
%                                 rtmodelLT(kk,imodelcomb,:,:,:) = reshape(thisRegGlm.predict(tbl(1:288,:)),[3,4,24]);
%                                 rtmodelTT(kk,imodelcomb,:) = thisRegGlm.predict(tbl(289:end,:));
%                                 lmRT{iexp}{isub, imodelcomb} =thisRegGlm;
%                                 BICRT(kk,imodelcomb) = thisRegGlm.ModelCriterion.BIC;
%                                 rtcoeff{imodelcomb}(kk,:) = thisRegGlm.Coefficients.Estimate;
%                                 rtT{imodelcomb}(kk,:) = thisRegGlm.Coefficients.tStat;
%                                 rtRMSE{imodelcomb}(kk,:) = nanstd(thisRegGlm.Residuals.LinearPredictor);

                                isLearnRel(imodelcomb) = isfield(modelsinfo{ilearnmodel},'contextual')&&modelsinfo{ilearnmodel}.contextual;
                                isDiffParams(imodelcomb) = (idiffpar==2 && iregcontextLT == iregcontextTT);
                                isDiffDQ(imodelcomb) = idiffBDQ==2;
                                isDiffB0(imodelcomb) = idiffB0 ==2;
                                whichConfModelLT(imodelcomb) = iregcontextLT;
                                whichConfModelTT(imodelcomb) = iregcontextTT;
                                whichLearnModel(imodelcomb) = ilearnmodel;
                            end
                        end
                    end
                end
            end
        end
    end
end


fitinfo.date = date;
save([resultsdir,'reg_conf_dual'],'glmCONF','BICCONF','AICCONF','confcoeff','conft','whichConfModelLT','whichConfModelTT','whichLearnModel','isDiffParams', 'isDiffB0','isDiffDQ','confmodelLT','confmodelTT','conf_mat','confRMSE','fitinfo','lmCONF')
