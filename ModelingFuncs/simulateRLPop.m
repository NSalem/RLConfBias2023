
function [out,out_sorted] = simulateRLPop(genparams,modelsinfo,nsess,ntrials,outcomemagnitudes)
%simulate RL behavior and hidden variables for a population given
%parameters and model
%genparams is cell array of size subs x learning models

n_sub = size(genparams,1);
nmodels = numel(modelsinfo);
out.pc = nan(n_sub,nmodels,nsess,ntrials*4);
out.dQ = nan(n_sub,nmodels,nsess,ntrials*4);
out.Qc = nan(n_sub,nmodels,nsess,ntrials*4);
out.Qu = nan(n_sub,nmodels,nsess,ntrials*4);
out.sigmaQ = nan(n_sub,nmodels,nsess,ntrials*4);
out.choice = nan(n_sub,nmodels,nsess,ntrials*4);
out.s=  nan(n_sub,nmodels,nsess,ntrials*4);
out.r=  nan(n_sub,nmodels,nsess,ntrials*4);
out.c=  nan(n_sub,nmodels,nsess,ntrials*4);
out.Q = nan(n_sub,nmodels,nsess,2,ntrials*4);
out.V = nan(n_sub,nmodels,nsess,ntrials*4);

out.choicePost = nan(n_sub,nmodels,112);
out.ss = nan(n_sub,nmodels,2,112);
out.dQ_post =  nan(n_sub,nmodels,112);
out.Q_c_post =  nan(n_sub,nmodels,112);
out.sigmaQ_post =  nan(n_sub,nmodels,112);
% 

out_sorted.pc = nan(n_sub,nmodels,nsess,4,ntrials);
out_sorted.choice = nan(n_sub,nmodels,nsess,4,ntrials);
out_sorted.s=  nan(n_sub,nmodels,nsess,4,ntrials);
out_sorted.r=  nan(n_sub,nmodels,nsess,4,ntrials);
out_sorted.c=  nan(n_sub,nmodels,nsess,4,ntrials);
out_sorted.Q = nan(n_sub,nmodels,nsess,2,4,ntrials);
out_sorted.V = nan(n_sub,nmodels,nsess,4,ntrials);
out_sorted.Qc = nan(n_sub,nmodels,nsess,4,ntrials);
out_sorted.Qu = nan(n_sub,nmodels,nsess,4,ntrials);

for isub = 1:size(genparams,1)
    for igenlearnmodel = 1:numel(modelsinfo)
        %%%%generate parameters
        paramstructgen = modelsinfo{igenlearnmodel};
        for iparam = 1:numel(modelsinfo{igenlearnmodel}.paramnames)
            thisParam = modelsinfo{igenlearnmodel}.paramnames{iparam};
            paramstructgen.(thisParam) = genparams{isub,igenlearnmodel}(iparam);
        end
        
        %%%%simulate RL for this subject and this model
        
        %%% generate outcome matrix and shuffle
        outcomes = Generate_Outcomes(24,outcomemagnitudes);
        condsorted = zeros(nsess,ntrials,4);
        for isess = 1:nsess
            condsorted(isess,:,:) = (ones(ntrials,1)*(1:4));
        end
        ss = repmat(nchoosek(1:8,2),4,1);
        permidcs = randperm(size(ss,1));
        ss = ss(permidcs,:);
        
        outcomes_all = zeros(nsess,2,ntrials*4);
        s = zeros(nsess,ntrials*4);
        for isess = 1:nsess
            s_sess = squeeze(condsorted(isess,:,:));
            permidcs = randperm(numel(s_sess));
            s_sess = s_sess(permidcs);
            s(isess,:) = s_sess;
            outcomes_sess = permute(outcomes,[1,3,2]);
            outcomes_sess = outcomes_sess(:,permidcs);
            outcomes_all(isess,:,:) = outcomes_sess;
        end
        
        %%% simulate single participant
        outSub = simulateRLSub(paramstructgen,s,outcomes_all,ss);
        
        %%%%fill matrices for all subjects and models
        out.pc(isub,igenlearnmodel,:,:) = outSub.pc;
        out.choice(isub,igenlearnmodel,:,:) = outSub.a;
        out.choicePost(isub,igenlearnmodel,:) = outSub.aa;
        out.s(isub,igenlearnmodel,:,:) = outSub.s;
        out.r(isub,igenlearnmodel,:,:) = outSub.r;
        out.c(isub,igenlearnmodel,:,:) = outSub.c;
        out.ss(isub,igenlearnmodel,:,:) = outSub.ss';
        out.Q(isub,igenlearnmodel,:,:,:)  = outSub.Q;
        out.dQ(isub,igenlearnmodel,:,:) = outSub.dQ;
        out.Qc(isub,igenlearnmodel,:,:) = outSub.Q_c;
        out.sigmaQ(isub,igenlearnmodel,:,:) = outSub.sigmaQ;
        out.V(isub,igenlearnmodel,:,:) =outSub.V;
        out.dQ_post(isub,igenlearnmodel,:) = outSub.dQ_post;
        out.Qc_post(isub,igenlearnmodel,:) = outSub.Q_c_post;
        out.Qu_post(isub,igenlearnmodel,:) = outSub.Q_uc_post;
        out.sigmaQ_post(isub,igenlearnmodel,:) = outSub.sigmaQ_post;
        out.V_post(isub,igenlearnmodel,:) = outSub.V_post;
        
        out.dQMat(isub,igenlearnmodel,:,:) = outSub.dQMat;
        out.QcMat(isub,igenlearnmodel,:,:) = outSub.QcMat;
        out.dQabsMat(isub,igenlearnmodel,:,:) = outSub.dQabsMat; 
        out.sigmaQMat(isub,igenlearnmodel,:,:) = outSub.sigmaQMat; 
        out.VMat(isub,igenlearnmodel,:,:) = outSub.VMat; 
        
        %%% sort per condition
        for isess = 1:nsess
            for icond = 1:4
                out_sorted.a(isub,igenlearnmodel,isess,icond,:) = outSub.a(isess,s(isess,:)==icond);
                out_sorted.pc(isub,igenlearnmodel,isess,icond,:) = outSub.pc(isess,s(isess,:)==icond);
                out_sorted.dQ(isub,igenlearnmodel,isess,icond,:) = outSub.dQ(isess,s(isess,:)==icond);
                out_sorted.Qc(isub,igenlearnmodel,isess,icond,:) = outSub.Q_c(isess,s(isess,:)==icond);
                out_sorted.Qu(isub,igenlearnmodel,isess,icond,:) = outSub.Q_uc(isess,s(isess,:)==icond);
                out_sorted.V(isub,igenlearnmodel,isess,icond,:) = outSub.V(isess,s(isess,:)==icond);
                out_sorted.sigmaQ(isub,igenlearnmodel,isess,icond,:) = outSub.sigmaQ(isess,s(isess,:)==icond);
            end
        end
        
        out_sorted.dQ_post(isub,igenlearnmodel,:) = outSub.dQ_post;
        out_sorted.Q_c_post(isub,igenlearnmodel,:) = outSub.Q_c_post;
        out_sorted.sigmaQ_post(isub,igenlearnmodel,:) = outSub.sigmaQ_post;
        out_sorted.V_post(isub,igenlearnmodel,:) = outSub.V_post;
        
    end
end

out.dQabs = abs(out.dQ);
out.dQabs_post = abs(out.dQ_post);

out_sorted.dQabs = abs(out_sorted.dQ);
out_sorted.dQabs_post = abs(out_sorted.dQ_post);

end