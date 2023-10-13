%% correlate Qc w/ all possible difficulty variables

load('Results\RLvars')

% dumQc = squeeze(nanmean(nanmean(nanmean(Q_c(11,:,:,:,:),3),4),5));

dumQc = squeeze(Q_c(11,:));
dumQGB =squeeze(Q(11,:,:,:,2,:)-Q(11,:,:,:,1,:))
dumQGB = dumQGB(:);
dumQCMU = squeeze(Q_c(11,:)-Q_uc(11,:));
dumdQAbs = abs(dumQCMU);
dumpc = squeeze(pc(11,:));


mats = [dumQCMU(:),dumdQAbs(:),dumQGB(:),dumpc(:)]';

figure()
for i = 1:size(mats,1)
subplot(2,2,i)
scatter(mats(i,:),dumQc);
r = corr(mats(i,:)',dumQc','rows','complete');

text(0,0,sprintf('r = %0.2f',r))
end