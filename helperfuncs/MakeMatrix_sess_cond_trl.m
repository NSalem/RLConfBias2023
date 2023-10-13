function data_out = MakeMatrix_sess_cond_trl(data,analyzeConfidence)
%%% Separates trials accross conditions for the learning task
%%% input:
%%%     data: (struct) Contains data structure of variables of size
%%%         subject x session x trial for learning task, and subject x trial
%%%         for the transfer task
%%%     analyzeConfidence: (bool) whether to include the confidence
%%%         variable
%%% output:
%%%     data_out: (struct) Contains data structure of size
%%%         subject x session x conditions x trial and subject x trial for the
%%%         transfer task

outMag = unique(abs(data.out));%outcome magnitudes
assert(numel(outMag)==2,'data contains more than two outcome magnitudes')

optP = [0.75,0.25,0.75,0.25,0.25,0.75,0.25,0.75]; %probability of high magnitude outcome per option
optU = [optP*outMag(2)+(1-optP)*outMag(1)]; %option utility (values in exp 1 are different but same order)
optU(5:8) = -optU(5:8);
subjects    = data.subjects;           % 3 sess
n_sub       = length(subjects);
n_sess      = size(data.con,2);
% n_trials    = size(data.con,3);
% expN        = data.expN;
for k_sub       = 1:n_sub
    for k_sess  = 1:n_sess
        con_sub = data.con(k_sub,k_sess,:);
        corr_sub = squeeze(data.cho(k_sub,k_sess,:)-1);
        out_sub = squeeze(data.out(k_sub,k_sess,:));
        out_good_sub = squeeze(data.out_good(k_sub,k_sess,:));
        cou_sub = squeeze(data.cou(k_sub,k_sess,:));
        cou_good_sub = squeeze(data.cou_good(k_sub,k_sess,:));
        RT_sub = squeeze(data.RT(k_sub,k_sess,:));
        RT_prev = [0;squeeze(data.RT(k_sub,k_sess,1:end-1))];
        
        if analyzeConfidence
            conf = data.conf(k_sub,k_sess,:);
            confN = [0;squeeze(data.conf(k_sub,k_sess,1:end-1))]; % confidence in previous trial REGARDLESS OF CONDITION
        end
        
        for k_cond = 1:4;
            cond = data.con(k_sub,k_sess,:);
            data_out.out_mat_sess(k_sub,k_sess,k_cond,:)    = out_sub(cond==k_cond);
            data_out.cou_mat_sess(k_sub,k_sess,k_cond,:)    = cou_sub(cond==k_cond);
            data_out.out_good_mat_sess(k_sub,k_sess,k_cond,:)    = out_good_sub(cond==k_cond);
            data_out.cou_good_mat_sess(k_sub,k_sess,k_cond,:)    = cou_good_sub(cond==k_cond);
            data_out.corr_mat_sess(k_sub,k_sess,k_cond,:)    = corr_sub(cond==k_cond);
            data_out.con_mat_sess(k_sub,k_sess,k_cond,:) = con_sub(cond==k_cond);
            data_out.RT_mat(k_sub,k_sess,k_cond,:)    = RT_sub(cond==k_cond);
            if analyzeConfidence
                data_out.conf_mat_sess(k_sub,k_sess,k_cond,:)    = conf(cond==k_cond);
                data_out.prevconf_mat_sess(k_sub,k_sess,k_cond,:)   = confN(cond==k_cond);
            end
        end
    end
    % recover post-learning prefs
    aa_sub = data.aa(k_sub,:);
    for k_symb = 1:8;
        S1 = data.ss(k_sub,:,1) == k_symb;
        S2 =  data.ss(k_sub,:,2)== k_symb;
        data_out.pref(k_sub,k_symb) = 100*(mean(data.aa(k_sub,S1)==1) + mean(data.aa(k_sub,S2)==2))./2;
%         pref2(k_sub,k_symb) =
%         100*(sum(data.aa(k_sub,S1)==1)+sum(data.aa(k_sub,S2)==2))/(sum(S1+S2));
%         % more similar to what makePostMats, but it's equivalent

        %%% get confidence when k_symb is present

        if analyzeConfidence
            data_out.confpostcon(k_sub,k_symb) = nanmean(data.conf_post(k_sub,(S1|S2)));
        end
                
        %%% get proportion of correct choice when k_symb is present
        optDiff = optU(data.ss(k_sub,:,2))-optU(data.ss(k_sub,:,1));
        whichCorrect = sign(optDiff);
        whichCorrect(whichCorrect==1) = 2;
        whichCorrect(whichCorrect==-1) = 1;
        
        correctChoice = aa_sub == whichCorrect;
        data_out.accPost(k_sub,k_symb) = sum((S1|S2) & correctChoice)./sum((S1|S2)& whichCorrect~=0);
        
    end
end
end