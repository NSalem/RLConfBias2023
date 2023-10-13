
load('Results\SimulationsSubjParams.mat')
% load('Results\SimulationsSubjParamsPerseveration.mat')

% rlVars = load('Results\RLvars_all.mat');
reg = load('Results\reg_conflogit_posttest_dqabs.mat');
% [prefMat_all,accMat_all,confdata,confMatPair_all] = makePostMats(rlVars.ss,rlVars.aa,rlVars.confpost)

idc = [6 5,2 1,8 7,4 3];
confmat = squeeze(nanmean(confChoiceMat(:,:,:,idc,idc)))*100;
confdata = reg.confmatChoice(:,idc,idc);
confdata(confdata==0) =NaN;
avgconf = squeeze(nanmean(confdata))'*100;

%%%  asymmetry in data

for isub = 1:size(confdata,1)
    confdatasub =squeeze(confdata(isub,:,:))'*100;
    
    for ir = 1:8
        for ic =1:8
            if isnan(confdatasub(ir,ic))
                confdatasub(ir,ic) = avgconf(ir,ic);
            end
        end
    end
    
    matdiffsub = (tril(confdatasub)-triu(confdatasub)');
    matdiffall(isub,:,:) = matdiffsub;
    matdiffsub(matdiffsub==0) = NaN;
    
    %     matdiffsub(isnan(matdiffsub)) = 0;
    asymsub = nanmean(matdiffsub(:));
    asymData(isub) = asymsub;
    
    %%% variance over chosen
    varconfc(isub) = var(nanmean(confdatasub,2));
    varconfu(isub) = var(nanmean(confdatasub,1));
    %
    y = confdatasub(:);
    fc = repmat(1:8,1,8);
    fu = repmat(1:8,8,1);
    fc = categorical(fc(:));
    fu = categorical(fu(:));
    tbl = table(confdatasub(:),fc(:),fu(:),'VariableNames',{'conf','fc','fu'});
    %
    %     anovan(y(~isnan(y)),{fc(~isnan(y)),fu(~isnan(y))})
    %     dum = fitlm(tbl);
    %     dum_fc = fitlm(tbl,'conf~1+fc');
    %     dum_fu = fitlm(tbl,'conf~1+fu');
end

confmat(confmat==0) =NaN;
%%% plot asymmetry of matrix per simulated model and data
asymModels = zeros(4,1);
for isub = 1:size(confmat,2)
    for imodel = 1:4
        mat = squeeze(confmat(imodel,isub,:,:))';
        
        for ir = 1:8
            for ic =1:8
                if isnan(mat(ir,ic))
                    mat(ir,ic) = squeeze(nanmean(confmat(imodel,:,ir,ic),2));
                end
            end
        end
        
        asymThisModel = (tril(mat)-triu(mat)');
        asymThisModel(asymThisModel==0) = NaN;

        matdiffallmodels(imodel,isub,:,:) = asymThisModel;
        
        %         asymThisModel(isnan(asymThisModel)) = 0;
        asymModels(imodel,isub) = nanmean(asymThisModel(:));
        
        %%% variance over chosen
        varconfcModels(imodel,isub) = var(nanmean(mat,2));
        varconfuModels(imodel,isub) = var(nanmean(mat,1));
    end
end

%% ratio of variance
figure()
varratio = (varconfc)./(varconfc+varconfu);
varratioModels = (varconfcModels)./(varconfcModels+varconfuModels);

m_data = nanmean(varratio);
se_data = nanstd(varratio)/sqrt(sum(~isnan(varratio)));

pirateplot(varratio,[0.5,0.5,0.5],0,1,10,'','','','')

m_sim = nanmean(varratioModels,2);
se_sim = nanstd(varratioModels,[],2)/sqrt(size(varratioModels,2));
hold on

colors = repmat([1,0.2,1],4,1).*(linspace(0.5,1,4)');
for imodel = 1:4
%     errorbar(imodel*0.5,m_sim(imodel),se_sim(imodel),'LineStyle','None','Marker','o')
    e(imodel) = errorbar(imodel*0.1+1-0.1*4/2-0.1/2,m_sim(imodel),se_sim(imodel),'LineStyle','None','Marker','.','Color',colors(imodel,:));
end

xlim([0,2])
% xticks(linspace(0,1,4))
% xticklabels({'Data','unbiased','Qc','\SigmaQ','V'})
leg = legend(e,{'-','Q_C','\Sigma Q','V'})
title(leg,'Simulated bias')
ylabel('confidence variance (C)/(C+U)')
