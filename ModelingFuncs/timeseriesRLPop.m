
function [out,out_sorted] = timeseriesRLPop(genparams,modelsinfo,s,a,r,c,ss,aa)
%genparams is cell array of size subs x learning models
nsub = size(genparams,1);
nsess = size(s,3);
ntrials = size(s,4)/4;

nmodels = numel(modelsinfo);
out.pc = nan(nsub,nmodels,nsess,ntrials*4);
out.dQ = nan(nsub,nmodels,nsess,ntrials*4);
out.Qc = nan(nsub,nmodels,nsess,ntrials*4);
out.sigmaQ = nan(nsub,nmodels,nsess,ntrials*4);
out.choice = nan(nsub,nmodels,nsess,ntrials*4);
out.s=  nan(nsub,nmodels,nsess,ntrials*4);
out.r=  nan(nsub,nmodels,nsess,ntrials*4);
out.c=  nan(nsub,nmodels,nsess,ntrials*4);
out.Q = nan(nsub,nmodels,nsess,2,ntrials*4);
out.V = nan(nsub,nmodels,nsess,ntrials*4);

out.choicePost = nan(nsub,nmodels,112);
out.ss = nan(nsub,nmodels,2,112);
out.dQ_post =  nan(nsub,nmodels,112);
out.Q_c_post =  nan(nsub,nmodels,112);
out.sigmaQ_post =  nan(nsub,nmodels,112);
% 

out_sorted.pc = nan(nsub,nmodels,nsess,4,ntrials);
out_sorted.choice = nan(nsub,nmodels,nsess,4,ntrials);
out_sorted.s=  nan(nsub,nmodels,nsess,4,ntrials);
out_sorted.r=  nan(nsub,nmodels,nsess,4,ntrials);
out_sorted.c=  nan(nsub,nmodels,nsess,4,ntrials);
out_sorted.Q = nan(nsub,nmodels,nsess,2,4,ntrials);
out_sorted.V = nan(nsub,nmodels,nsess,4,ntrials);

for isub = 1:size(genparams,1)
    for igenlearnmodel = 1:numel(modelsinfo)
        %%%%generate parameters
        paramstructgen = modelsinfo{igenlearnmodel};
        for iparam = 1:numel(modelsinfo{igenlearnmodel}.paramnames)
            thisParam = modelsinfo{igenlearnmodel}.paramnames{iparam};
            paramstructgen.(thisParam) = genparams{isub,igenlearnmodel}(iparam);
        end
        
        %%%%simulate RL for this subject and this model
        
        s_sub = squeeze(s(isub,igenlearnmodel,:,:));
        a_sub = squeeze(a(isub,igenlearnmodel,:,:));
        r_sub = squeeze(r(isub,igenlearnmodel,:,:));
        c_sub = squeeze(c(isub,igenlearnmodel,:,:));
        ss_sub = squeeze(ss(isub,igenlearnmodel,:,:));
        aa_sub = aa(isub,igenlearnmodel,:);
        
        %%% simulate single participant
        outSub = timeSeriesRLSub(paramstructgen,s_sub,a_sub,r_sub,c_sub,ss_sub',aa_sub);
        
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
        out.Qu(isub,igenlearnmodel,:,:) = outSub.Q_uc;
        out.sigmaQ(isub,igenlearnmodel,:,:) = outSub.sigmaQ;
        out.V(isub,igenlearnmodel,:,:) =outSub.V;
        out.dQ_post(isub,igenlearnmodel,:) = outSub.dQ_post;
        out.Q_post(isub,igenlearnmodel,:,:) = outSub.Q_post;
        out.Qc_post(isub,igenlearnmodel,:) = outSub.Q_c_post;
        out.Qu_post(isub,igenlearnmodel,:) = outSub.Q_uc_post;
        out.sigmaQ_post(isub,igenlearnmodel,:) = outSub.sigmaQ_post;
        out.V_post(isub,igenlearnmodel,:) = outSub.V_post;
        out.pc_post(isub,igenlearnmodel,:) = outSub.pc_post;
        
        %%% sort per condition
        for isess = 1:nsess
            for icond = 1:4
                out_sorted.a(isub,igenlearnmodel,isess,icond,:) = outSub.a(isess,s_sub(isess,:)==icond);
                out_sorted.pc(isub,igenlearnmodel,isess,icond,:) = outSub.pc(isess,s_sub(isess,:)==icond);
                out_sorted.dQ(isub,igenlearnmodel,isess,icond,:) = outSub.dQ(isess,s_sub(isess,:)==icond);
                out_sorted.Qc(isub,igenlearnmodel,isess,icond,:) = outSub.Q_c(isess,s_sub(isess,:)==icond);
                out_sorted.Qu(isub,igenlearnmodel,isess,icond,:) = outSub.Q_uc(isess,s_sub(isess,:)==icond);
                out_sorted.V(isub,igenlearnmodel,isess,icond,:) = outSub.V(isess,s_sub(isess,:)==icond);
                out_sorted.sigmaQ(isub,igenlearnmodel,isess,icond,:) = outSub.sigmaQ(isess,s_sub(isess,:)==icond);
            end
        end
        
        out_sorted.dQ_post(isub,igenlearnmodel,:) = outSub.dQ_post;
        out_sorted.Q_c_post(isub,igenlearnmodel,:) = outSub.Q_c_post;
        out_sorted.sigmaQ_post(isub,igenlearnmodel,:) = outSub.sigmaQ_post;
        out_sorted.V_post(isub,igenlearnmodel,:) = outSub.V_post;
                
        out.dQMat(isub,igenlearnmodel,:,:) = outSub.dQMat;
        out.QcMat(isub,igenlearnmodel,:,:) = outSub.QcMat;
        out.dQabsMat(isub,igenlearnmodel,:,:) = outSub.dQabsMat; 
        out.sigmaQMat(isub,igenlearnmodel,:,:) = outSub.sigmaQMat; 
        out.VMat(isub,igenlearnmodel,:,:) = outSub.VMat; 
        
    end
end
out.dQabs = abs(out.dQ);
out.dQabs_post = abs(out.dQ_post);
end