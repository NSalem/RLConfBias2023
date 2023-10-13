function [out] = timeseriesRLSub(paramstructgen,s,a,r,c,ss,aa)
%%% simulate hidden variables given RL model and actions

nsess = size(s,1);
ntrials = size(s,2);

out.Q = nan(nsess,2,ntrials);
out.V = nan(nsess,ntrials);
out.a = a;
out.r = nan(nsess,ntrials);
out.c = nan(nsess,ntrials);
out.s = s;
out.pc= nan(nsess,ntrials);
out.conf = nan(nsess,ntrials);

out.ss = ss;

for isess = 1:nsess
    s_sess = squeeze(s(isess,:));
    a_sess = squeeze(a(isess,:));
    r_sess = squeeze(r(isess,:));
    c_sess = squeeze(c(isess,:));
    
    if isess == nsess
        dum_ss = out.ss;
    else
        dum_ss = [];
    end
    
    [Q_sess,V_sess,pc_sess,PE_sess,ppc,pll,dQ_post,V_post,Q_post,Qf,Vf] = Computational_TimeSeries_QLearner_Unsorted(paramstructgen,s_sess(:),a_sess(:),r_sess(:),c_sess(:),aa,ss);
    out.Q(isess,:,:) = Q_sess(:,1:ntrials);
    out.V(isess,:) = V_sess(1:ntrials);
    out.a(isess,:,:) = a_sess;
    out.r(isess,:) = r_sess;
    out.c(isess,:) = c_sess;
    out.pc(isess,:) = pc_sess;
    out.aa = aa;
    %conditions in order of presentation
end

clear condsorted
out.Q_post = Q_post;
out.V_post = V_post;
out.dQ = squeeze(out.Q(:,2,:)-out.Q(:,1,:));
out.Q_c = nan(nsess,ntrials);
out.Q_uc = nan(nsess,ntrials);
out.dQ_post = dQ_post;
out.Q_c_post = nan(1,length(Q_post));
out.Q_uc_post = nan(1,length(Q_post));
out.pc_post = nan(1,length(Q_post));

for isess = 1:nsess
    for itrl = 1:ntrials
        if isnan(out.a(isess,itrl))
            out.Q_c(isess,itrl) = NaN;
            out.Q_uc(isess,itrl) = NaN;
        else
            out.Q_c(isess,itrl) = out.Q(isess,out.a(isess,itrl),itrl);
            out.Q_uc(isess,itrl) = out.Q(isess,3-out.a(isess,itrl),itrl);
        end
    end
end
out.sigmaQ = out.Q_c+out.Q_uc;

optP = [0.75,0.25,0.75,0.25,0.25,0.75,0.25,0.75];
optU = [optP*1+(1-optP)*0.1]; %option utility (values in exp 1 are different but same order)
optU(5:8) = -optU(5:8);

for itrl = 1:length(Q_post)
    out.Q_c_post(:,itrl) = out.Q_post(out.aa(itrl),itrl);
    out.Q_uc_post(:,itrl) = out.Q_post(3-out.aa(itrl),itrl);
    
    if diff(optU(ss(itrl,:)))>0
        ig = 2;
    elseif  diff(optU(ss(itrl,:)))<0
        ig = 1;
    else
        ig = [];
    end
    if ~isempty(ig)
        
        if ig == aa(itrl)
            out.pc_post(itrl) = pll(itrl);
        else
            out.pc_post(itrl) = 1-pll(itrl);
        end
    end
end
out.sigmaQ_post = out.Q_c_post+out.Q_uc_post;

for k_symb = 1:8; %chosen
    for k_symb2 = 1:8 %unchosen
        out.dQMat(k_symb2,k_symb) = Qf(k_symb2)-Qf(k_symb);
        out.dQabsMat(k_symb2,k_symb) = abs(Qf(k_symb2)-Qf(k_symb));
        out.QcMat(k_symb2,k_symb) = Qf(k_symb2);
        out.sigmaQMat(k_symb2,k_symb) =  Qf(k_symb2)+Qf(k_symb);
        out.VMat(k_symb2,k_symb) = nanmean(Vf([k_symb2,k_symb]));
    end
end
end
