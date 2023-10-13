function [prefMat,accMat,confChoiceMat,confChoicePairMat,pref,accPost,confPost]= makePostMats(ss,aa,conf)
%%%ss: subj x trial x option  
%%%aa: subj x trial
%%%conf: subj x trial
optP = [0.75,0.25,0.75,0.25,0.25,0.75,0.25,0.75];
optU = [optP*1+(1-optP)*0.1]; %option utility (values in exp 1 are different but same order)
optU(5:8) = -optU(5:8);            
accMat = NaN(size(ss,1),8,8);
prefMat = NaN(size(ss,1),8,8);
confChoiceMat = NaN(size(ss,1),8,8); 
confChoicePairMat= NaN(size(ss,1),8,8);
confPost = NaN(size(ss,1),8);
pref = NaN(size(ss,1),8);
for isub = 1:size(ss,1)
    ss_sub = squeeze(ss(isub,:,:));
    aa_sub = squeeze(aa(isub,:));
    for k_symb = 1:8; 
        SC1 = ss_sub(:,1) == k_symb;
        SC2 = ss_sub(:,2)== k_symb;
%         prefOld(isub,k_symb) = 100*nanmean([nanmean(aa_sub(SC1)==1), nanmean(aa_sub(SC2)==2)]);
        pref(isub,k_symb) = 100*(sum(aa_sub(SC1)==1)+sum(aa_sub(SC2)==2))/(sum(SC1+SC2));
%         pref3(isub,k_symb) = 100*(mean(aa(isub,SC1)==1) + mean(aa(isub,SC2)==2))./2;
        optDiff = optU(ss_sub(:,2))-optU(ss_sub(:,1));
        whichCorrect = sign(optDiff);
        whichCorrect(whichCorrect==1) = 2;
        whichCorrect(whichCorrect==-1) = 1;
%         
        correctChoice = aa_sub == whichCorrect;
        accPost(isub,k_symb) = sum((SC1|SC2) & correctChoice' & whichCorrect'~=0)./sum((SC1|SC2)& whichCorrect'~=0);

        if ~isempty(conf)
            confPost(isub,k_symb) = nanmean(conf(isub,find([SC1|SC2])));
        end
        for k_symb2 = 1:8 
            
            S1 = (ss_sub(:,1) == k_symb) & (ss_sub(:,2) == k_symb2); %where pair of symbols appear in order 1
            S2 = (ss_sub(:,2) == k_symb) & (ss_sub(:,1) == k_symb2); %where pair of symbols appear in order 2
            nChose1 = nansum(aa_sub(S1)==1)+nansum(aa_sub(S2)==2); %chosen
            optDiff = optU(k_symb)-optU(k_symb2);
            if ~isempty(conf)
                confChoiceMat(isub,k_symb2,k_symb) = (nanmean(conf(isub,(S1& aa_sub' ==1 )|(S2&aa_sub' ==2))));            
                confChoicePairMat(isub,k_symb2,k_symb) = (nanmean(conf(isub,S1|S2)));
            end
            if optDiff~=0
                accMat(isub,k_symb,k_symb2) = nChose1.*(optDiff>0)+...
                    (4-nChose1).*(optDiff<0);
            else
                
            end
            
            if k_symb~=k_symb2
                prefMat(isub,k_symb,k_symb2) = nChose1;
            end
        end
    end
end
end