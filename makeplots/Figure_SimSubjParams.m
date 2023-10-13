addpath('helperfuncs')
load('Results\SimulationsSubjParams.mat')

col_mat = [0,1,0;...
    0,1,1;...
    1,0,0;...
    1,0,1];


Xfill = [1:24,fliplr(1:24)];


% close()
newIdc = [6,5,2,1,8,7,4,3];
prefMat= prefMat(:,:,newIdc,newIdc);
accMat=  accMat(:,:,newIdc,newIdc);
confChoiceMat = confChoiceMat(:,:,:,newIdc,newIdc);
confChoicePairMat = confChoicePairMat(:,:,:,newIdc,newIdc);

prefMat= squeeze(nanmean(prefMat));
accMat = squeeze(nanmean(accMat));
confChoiceMat = squeeze(nanmean(confChoiceMat));
confChoicePairMat = squeeze(nanmean(confChoiceMat));

[h1,h2]= makePlotsSims(outSorted,outConfLSorted,prefMat,accMat,confChoiceMat)

saveas(h1,['Plots/simSubjParams/simBehavLT.svg'])
saveas(h2,['Plots/simSubjParams/simBehavTT.svg'])
% close()
%%
% thispc =  mean(squeeze(nanmean(nanmean(nanmean(pc_sorted,4),5),3)));
thispc = squeeze(nanmean(nanmean(outSorted.pc,3),5));
for imodel = 1:numel(idmodels)
	meanConfCond = squeeze(nanmean(nanmean(nanmean(outConfLSorted(:,:,imodel,:,:,:),1),4),6));
    confBiasLT{imodel} = nanmean(squeeze(mean(meanConfCond(:,[1,2])-meanConfCond(:,[3,4]),2)),2);
    overconfLT{imodel} = nanmean(meanConfCond,2)-nanmean(thispc,2).*100;
end

newIdc = [6,5,2,1,8,7,4,3]; %rearange order of post-test stimuli
% accTT = squeeze(nanmean(accPost));
accTT = accPost(:,newIdc).*100;

for imodel= 1:numel(idcTT)
    confTT = conf_matTT{imodel};
    confTT = confTT(:,newIdc);
    confTTAll{imodel} = confTT;
    overconfTT{imodel} = nanmean(confTT-accTT,2);
    confBiasTT{imodel} = nanmean(confTT(:,[3,4,7,8])-confTT(:,[1,2,5,6]),2);
end

%%% plot asymmetry of matrix per model
asymAll = zeros(4,1);
for imodel = 1:4
    mat = squeeze(nanmean(confChoiceMat(imodel,:,:,:),2))'*100;
    asymThisModel = abs(tril(mat)-tril(mat'));
    asymThisModel(asymThisModel==0) = NaN;
    asymAll(imodel) = nanmean(asymThisModel(:));
end