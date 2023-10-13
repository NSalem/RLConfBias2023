function [out] = simulateRLSub(paramstructgen,s,OUT,ss)
%%% simulate model actions and hidden variables for single participant
%%% input:
%     paramstructgen: struct of parameters and specifications for the model
%     s: array of states
%     OUT: array of outcomes allocated on each side
%     ss: array of transfer task symbols
%%% output:
%     out: struct with fields
%         a: array of actions in learning task
%         r: array of obtained outcomes in learning task
%         c: array of forgone outcomes in learning task 
%         pc: array of probability of correct choice in learning task
%         Q: array of Q values across learning task
%         V: array of context values accross learning task
%         conf: placeholder array (empty) for confidence across learning task
    
nsess = size(OUT,1);
ntrials = size(OUT,3);

out.Q = nan(nsess,2,ntrials);
out.V = nan(nsess,ntrials);
out.a = nan(nsess,ntrials);
out.r = nan(nsess,ntrials);
out.c = nan(nsess,ntrials);
out.s = nan(nsess,ntrials);
out.pc= nan(nsess,ntrials);
out.conf = nan(nsess,ntrials);

out.ss = ss;

for isess = 1:nsess
    s_sess = squeeze(s(isess,:));
    OUTsess = squeeze(OUT(isess,:,:));
    COMP = (1-mod(s_sess,2));
    
    if isess == nsess
        dum_ss = out.ss;
    else
        dum_ss = [];
    end
    [Q_sess,V_sess,pc_sess,a_sess,r_sess,c_sess,aa,dQ_post,V_post,Q_post,Qf,Vf] = Computational_Simus_QLearner_Unsorted(paramstructgen,s_sess,dum_ss,OUTsess);
    
    out.Q(isess,:,:) = Q_sess(:,1:ntrials);
    out.V(isess,:) = V_sess(1:ntrials);
    out.a(isess,:,:) = a_sess;
    out.r(isess,:) = r_sess;
    out.c(isess,:) = c_sess;
    out.s(isess,:) = s_sess;
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
for itrl = 1:length(Q_post)
    out.Q_c_post(:,itrl) = out.Q_post(out.aa(itrl),itrl);
    out.Q_uc_post(:,itrl) = out.Q_post(3-out.aa(itrl),itrl);
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
