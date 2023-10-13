function saveModelVars(dataFile,fitFiles,outFile)
%%% generate hidden variables from behavior and  fitted parameters and saves
%%% results
%%% dataFile: name of data file containing behavior data
%%% fitFiles: file name prefix of files containing the fit data (minus the
%%% number indicating the experiment)
%%% outFile: file name for saving data
%%%
%%%    Q: time course of Q-values for all models and subjects
%%%     in a matrix of size model x subject x session x context x action x trial
%%%    V contains the context values (where applicable) in a matrix of size
%%%    model x subject x session x context x trial
%%%    Q_post_symb: the final Q-values for each symbol (model x
%%%    subject x symbol)
%%%    pc contains the time course of probability of correct choice
%%%    (model x subject x session x context x trial)
%%%    Q_c: time course of Q-value of the option that was chosen (model x subject x session x context x trial)
%%%    Q_u: time course of Q-value of the option that was NOT chosen (model x subject x session x context x trial)

% assert(numel(dataFile)==numel(fitFiles));

%%% load data %%%
exps = load([dataFile]);
exps_conf = exps.data_all;

%%% make lis tof number of sessions and trials per condition in each experiment %%%%
load([fitFiles,'1','.mat'],'modelsinfo');
models = 1:numel(modelsinfo);
ntrialsExp = zeros(size(exps_conf));
nsessExp = zeros(size(exps_conf));
ntrialspostExp = zeros(size(exps_conf));
nsub = 0;
for iexp = 1:numel(exps_conf)
    thisExp = MakeMatrix_sess_cond_trl(exps_conf(iexp),0);
    ntrialsExp(iexp) = size(thisExp.con_mat_sess,4);
    nsessExp(iexp) = size(thisExp.con_mat_sess,2);
    nsub=nsub+size(thisExp.corr_mat_sess,1);
    ntrialspostExp(iexp) = size(exps_conf(iexp).ss,2)    
end

maxTrialsExp = max(ntrialsExp);
maxSessExp = max(nsessExp);
maxTrialsPost = max(ntrialspostExp);


%%% model x sub x sess x cond x action x trial
Q = nan(numel(models),nsub,maxSessExp,4,2,maxTrialsExp);
%%% model x sub x sess x cond x trial
correct =  nan(nsub,maxSessExp,4,maxTrialsExp);
conf =  nan(nsub,maxSessExp,4,maxTrialsExp);
confprev =  nan(nsub,maxSessExp,4,maxTrialsExp);

dQ = nan(numel(models),nsub,maxSessExp,4,maxTrialsExp);
V = nan(numel(models),nsub,maxSessExp,4,maxTrialsExp);
Q_c = nan(numel(models),nsub,maxSessExp,4,maxTrialsExp);
Q_uc = nan(numel(models),nsub,maxSessExp,4,maxTrialsExp);
pc =  nan(numel(models),nsub,maxSessExp,4,maxTrialsExp);
pa =  nan(numel(models),nsub,maxSessExp,4,maxTrialsExp);

%%% model x sub x trials (transfer task)
ppc = nan(numel(models),nsub,8);
pll = nan(numel(models),nsub,maxTrialsPost);
Q_post_symb = nan(numel(models),nsub,8);
dQ_post = nan(numel(models),nsub,maxTrialsPost);
Q_c_post = nan(numel(models),nsub,maxTrialsPost);
Q_uc_post = nan(numel(models),nsub,maxTrialsPost);
QG_post = nan(numel(models),nsub,maxTrialsPost);
QB_post = nan(numel(models),nsub,maxTrialsPost);

confpost = nan(nsub,maxTrialsPost);

%%% model x sub x symb x symb 
dQMat = nan(numel(models),nsub,8,8);
absDQMat = nan(numel(models),nsub,8,8);
sigmaQMat = nan(numel(models),nsub,8,8);


kk = 0;
sub_exp = []; %Experiment number for each subject
sub_exp_name = {}; %Experiment folder name for each subject
sub_ntrials = [];
sub_nsess = [];

for iexp = 1:numel(exps_conf)
    infostr.subid{iexp} = zeros(size(exps_conf(iexp).subjects)); %list of subjects per exp
    %%% load model fits for experiment
    load([fitFiles,num2str(iexp),'.mat'],'parametersLPP','modelsinfo','LAME','ll');
    models = 1:size(parametersLPP,2);
    
    %%% loop over subjects to get hidden variables
    
    if isfield(exps_conf(iexp),'conf') && ~isempty(exps_conf(iexp).conf)
        analyzeConf = 1;
    else 
        analyzeConf = 0;
    end
    thisExp = MakeMatrix_sess_cond_trl(exps_conf(iexp),analyzeConf);
    
    for isub = 1:size(exps_conf(iexp).corr,1)
        
        %%% preallocate variables for participant %%%
        %%% with nans size maxSessions x conditions x maxTrials %%%
        dQ_sub = nan(maxSessExp,4,maxTrialsExp);
        V_sub = nan(maxSessExp,4,maxTrialsExp);
        pc_sub = nan(maxSessExp,4,maxTrialsExp);
        pa_sub = nan(maxSessExp,4,maxTrialsExp);
        Q_sub = nan(maxSessExp,4,2,maxTrialsExp);
        correct_sub = nan(maxSessExp,4,maxTrialsExp);
        conf_sub = nan(maxSessExp,4,maxTrialsExp);
        confprev_sub = nan(maxSessExp,4,maxTrialsExp);

        kk = kk+1;
        
        ss_sub = squeeze(exps_conf(iexp).ss(isub,:,:));
        aa_sub = squeeze(exps_conf(iexp).aa(isub,:));
        
        
        LAME_allexp(kk,:)=  LAME(isub,:);
        ll_allexp(kk,:)=  ll(isub,:);
        
        ntrials = size(thisExp.corr_mat_sess,4);
        nsess = size(exps_conf(iexp).corr,2);
        imodelcomb = 0;
        for ilearnmodel = models
            params_sub = parametersLPP{isub,ilearnmodel};
            params_allexp{kk,ilearnmodel} =  parametersLPP{isub,ilearnmodel};
            for isess = 1:nsess
                con_sub_sess = squeeze(thisExp.con_mat_sess(isub,isess,:,:))';
                %                 cho_sub_sess = nan(size(4,maxTrialsExp));
                cho_sub_sess = squeeze(thisExp.corr_mat_sess(isub,isess,:,:))'+1;
                out_sub_sess = squeeze(thisExp.out_mat_sess(isub,isess,:,:))';
                cou_sub_sess = squeeze(thisExp.cou_mat_sess(isub,isess,:,:))';
                
                if (analyzeConf)
                    conf_sub_sess =  squeeze(thisExp.conf_mat_sess(isub,isess,:,:))';
                    confprev_sub_sess =  squeeze(thisExp.prevconf_mat_sess(isub,isess,:,:))';
                end
                
                paramstruct = modelsinfo{ilearnmodel};
                for iparam = 1:numel(modelsinfo{ilearnmodel}.paramnames)
                    paramstruct.(modelsinfo{ilearnmodel}.paramnames{iparam}) = params_sub(iparam);
                end
                %                 if size(cho_sub_sess,1)<size(cho_sub_sess,2)
                %                     cho_sub_sess = cho_sub_sess';
                %                 end
                correct_sub(isess,:,1:ntrials) =cho_sub_sess'-1;
                
                if analyzeConf
                    conf_sub(isess,:,1:ntrials) =conf_sub_sess';
                    confprev_sub(isess,:,1:ntrials) =confprev_sub_sess';
                end
                %                 correct(kk,isess,:,1:ntrials) = cho_sub_sess'-1;
                
                %% Calculate hidden variables
                [Q_sub_sess,V_sub_sess,pc_sub_sess,PE,ppc_sub,pll_sub,dQ_post_sub,V_post_sub,Q_post,Qf,Vf]  = Computational_TimeSeries_QLearner(paramstruct,con_sub_sess,cho_sub_sess(:),out_sub_sess(:),cou_sub_sess(:),aa_sub,ss_sub);
                dQ_sub_sess = (squeeze(Q_sub_sess(:,2,1:ntrials)-Q_sub_sess(:,1,1:ntrials)));
                pa_sub_sess = pc_sub_sess;
                pa_sub_sess(cho_sub_sess' == 1) = 1-pa_sub_sess(cho_sub_sess' == 1);
                
                
%                 dQ_sub(isess,:,1:ntrials) = dQ_sub_sess;
                V_sub(isess,:,1:ntrials) = V_sub_sess(:,1:ntrials);
                pc_sub(isess,:,1:ntrials) = pc_sub_sess;
                pa_sub(isess,:,1:ntrials)= pa_sub_sess;
                Q_sub(isess,:,:,1:ntrials) = Q_sub_sess(:,:,1:ntrials);
                
                
                %                 dQ(ilearnmodel,kk,isess,1:4,1:ntrials) = dQ_sub_sess;
                %                 V(ilearnmodel,kk,isess,1:4,1:ntrials) = V_sub_sess(:,1:ntrials);
                %                 pc(ilearnmodel,kk,isess,:,:) = pc_sub_sess;
                %                 pa(ilearnmodel,kk,isess,:,:) = pa_sub_sess;
                %                 Q(ilearnmodel,kk,isess,:,:,:) = Q_sub_sess(:,:,1:ntrials);
                
                for icon = 1:4
                    for itrl = 1:ntrials
                        if isnan(cho_sub_sess(itrl,icon))
                            Q_c(ilearnmodel,kk,isess,icon,itrl) = NaN;
                            Q_uc(ilearnmodel,kk,isess,icon,itrl) = NaN;
                        else
                            Q_c(ilearnmodel,kk,isess,icon,itrl) = Q_sub_sess(icon,cho_sub_sess(itrl,icon),itrl);
                            Q_uc(ilearnmodel,kk,isess,icon,itrl) = Q_sub_sess(icon,3-cho_sub_sess(itrl,icon),itrl);
                        end
                    end
                end
            end
%             dQ(ilearnmodel,kk,:,:,:) = dQ_sub;
            V(ilearnmodel,kk,:,:,:) = V_sub;
            pc(ilearnmodel,kk,:,:,:) = pc_sub;
            pa(ilearnmodel,kk,:,:,:) = pa_sub;
            Q(ilearnmodel,kk,:,:,:,:) = Q_sub;
            ppc(ilearnmodel,kk,:) = ppc_sub;
            pll(ilearnmodel,kk,:) = pll_sub;
            correct(kk,:,:,:) = correct_sub;
            ss(kk,:,:) = ss_sub;
            aa(kk,:) = aa_sub;
            pref(kk,:) = thisExp.pref(isub,:);
            %%% transfer task variables
            Q_post_symb(ilearnmodel,kk,:) = Qf;
            
            
            if analyzeConf
                conf(kk,:,:,:,:) = conf_sub;
                confprev(kk,:,:,:,:) = confprev_sub;
                confpost(kk,:) = exps_conf(iexp).conf_post(isub,:);
            else 
                conf(kk,:,:,:,:) = nan.*correct_sub;
                confprev(kk,:,:,:,:) = nan.*correct_sub;
                confpost(kk,:) = nan.*aa_sub;
            end
            
            %%% get real utility
            optP = [0.75,0.25,0.75,0.25,0.25,0.75,0.25,0.75];
            optU = [optP*1+(1-optP)*0.1]; %option utility (values in exp 1 are different but same order)
            optU(5:8) = -optU(5:8);
            %             optU = optU([8,7,2,1,6,5,4,3]);
            
            %             id = [3,7,1,5,4,8,2,6];
            
            dQ_post(ilearnmodel,kk,:) = dQ_post_sub;
            V_post(ilearnmodel,kk,:) = V_post_sub;
            for itrl = 1:numel(aa_sub);
                if isnan(aa_sub(itrl))
                    Q_c_post(ilearnmodel,kk,itrl) = NaN;
                    Q_uc_post(ilearnmodel,kk,itrl) = NaN;
                    QG_post(ilearnmodel,kk,itrl) = NaN;
                    QB_post(ilearnmodel,kk,itrl) = NaN;
                else
                    Q_c_post(ilearnmodel,kk,itrl) = Q_post(aa_sub(itrl),itrl);
                    Q_uc_post(ilearnmodel,kk,itrl) = Q_post(3-aa_sub(itrl),itrl);
                    
                    if diff(optU(ss_sub(itrl,:)))>0
                        ig = 2;
                    elseif  diff(optU(ss_sub(itrl,:)))<0
                        ig = 1;
                    else
                        ig = [];
                    end
                    if ~isempty(ig)
                        QG_post(ilearnmodel,kk,itrl) = Q_post(ig,itrl);
                        QB_post(ilearnmodel,kk,itrl) = Q_post(3-ig,itrl);
                        
                        if ig == aa_sub(itrl)
                            pc_post(ilearnmodel,kk,itrl) = pll(ilearnmodel,kk,itrl);
                        else
                            pc_post(ilearnmodel,kk,itrl) = 1-pll(ilearnmodel,kk,itrl);
                        end
                        
                    else
                        pc_post(ilearnmodel,kk,itrl) = NaN;
                        QG_post(ilearnmodel,kk,itrl) = NaN;
                        QB_post(ilearnmodel,kk,itrl) = NaN;
                    end
                end
            end
            sigmaQ_post_sub = sum(Q_post);
            
            %         pair_mat = [6 5;2 1;8 7;4 3;];
            for k_symb = 1:8; %chosen
                for k_symb2 = 1:8 %unchosen
                    
                    dQMat(ilearnmodel,kk,k_symb2,k_symb) = Qf(k_symb2)-Qf(k_symb);
                    absDQMat(ilearnmodel,kk,k_symb2,k_symb) = abs(Qf(k_symb2)-Qf(k_symb));
                    sigmaQMat(ilearnmodel,kk,k_symb2,k_symb) =  Qf(k_symb2)+Qf(k_symb);
                    VMat(ilearnmodel,kk,k_symb2,k_symb) = nanmean(Vf([k_symb2,k_symb]));
                end
            end
            
%             clear conf_sub confprev_sub
        end
        sub_exp = [sub_exp,iexp];
        sub_exp_name = [sub_exp_name,exps_conf(iexp).exp];
        
        sub_ntrials =[sub_ntrials,ntrials];
        sub_nsess = [sub_nsess,nsess];
    end
    clear parametersLPP
    infostr.modelsinfo{iexp} = modelsinfo;
    infostr.subid{iexp} = [infostr.subid{iexp},kk];
    infostr.subexp = sub_exp;
    infostr.subexpname =  sub_exp_name;
    infostr.sub_ntrials = sub_ntrials;
    infostr.sub_nsess = sub_nsess;
end

dQ = Q_c - Q_uc;
dQabs = abs(Q_c - Q_uc);

% correct = permute(correct,[1,2,4,3]);
infostr.date = date;
infostr.models = models;
save(['Results\',outFile],'Q','dQ','dQabs','V','pc','pa','Q_c','Q_uc','ss','aa','pc_post','dQ_post','Q_c_post','Q_uc_post','V_post','correct','conf','confprev','ppc','pll','infostr','params_allexp','modelsinfo','dQMat','sigmaQMat','absDQMat','VMat','Q_post_symb','QG_post','QB_post','pc_post','LAME_allexp','ll_allexp','confpost','pref');

end