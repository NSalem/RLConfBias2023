clear all 
close all force

% cd('..')
load('Results\RLvars_all')
reg = load('Results\reg_conflogit_posttest_dqabs.mat');
ilearnmodel = 11;

% load('Results\RLvars_ExtraContext_all.mat');
% reg = load('Results\reg_conflogitPerseveration_posttest_dqabs.mat');
% ilearnmodel = 3;


idc = [6 5,2 1,8 7,4 3];

lg={'\color{green}P G_{75}','\color{green}P G_{25}','\color{cyan}C G_{75}','\color{cyan}C G_{25}',...
    '\color{red}P L_{25}','\color{red}P L_{75}','\color{magenta}C L_{25}','\color{magenta}C L_{75}'};

lg={lg{idc}};

reg.confmodelmatChoice = reg.confmodelmatChoice(:,:,idc,idc);


figure()
subplot(3,3,1)

myMatMean = squeeze(nanmean(absDQMat(ilearnmodel,:,idc,idc),2));
myMatMean(myMatMean==0) = NaN;
b = imagesc(myMatMean);
set(b,'AlphaData',~isnan(myMatMean));
    
set(gca,'FontName','Arial','FontSize',12)
xlabel('Unchosen')
ylabel('Chosen')

% colormap(jet)
h = colorbar();
title(h,'|\DeltaQ|')
xticklabels(lg);yticklabels(lg);xtickangle(90)
xticks(1:8);yticks(1:8)

subplot(3,3,2)

myMatMean = (squeeze(nanmean(Q_post_symb(ilearnmodel,:,idc),2)).*ones(1,8));
myMatMean(myMatMean==0) = NaN;
b = imagesc(myMatMean);
set(b,'AlphaData',~isnan(myMatMean));
set(b,'AlphaData',~eye(size(myMatMean)));

yticks(1:8);yticklabels(lg)
% xticks(lg);xtickangle(90)
h = colorbar()
title(h,'Q_c')
xticklabels(lg);yticklabels(lg);xtickangle(90)
xticks(1:8);yticks(1:8)
set(gca,'FontName','Arial','FontSize',12)
xlabel('Unchosen')
ylabel('Chosen')

subplot(3,3,5)

myMatMean = squeeze(nanmean(sigmaQMat(ilearnmodel,:,idc,idc),2));
myMatMean(myMatMean==0) = NaN;
b = imagesc(myMatMean);
set(b,'AlphaData',~isnan(myMatMean));
set(b,'AlphaData',~eye(size(myMatMean)));


% colormap(jet)
h = colorbar();
title(h,'\SigmaQ')
xticklabels(lg);yticklabels(lg)
xticks(1:8);yticks(1:8);xtickangle(90)
set(gca,'FontName','Arial','FontSize',12)
xlabel('Unchosen')
ylabel('Chosen')

subplot(3,3,8)

myMatMean = squeeze(nanmean(VMat(ilearnmodel,:,idc,idc),2));
myMatMean(myMatMean==0) = NaN;
b = imagesc(myMatMean);
set(b,'AlphaData',~isnan(myMatMean));
set(b,'AlphaData',~eye(size(myMatMean)));

% colormap(jet)
h = colorbar();
title(h,'V')
xticklabels(lg);yticklabels(lg)
xticks(1:8);yticks(1:8);xtickangle(90)
set(gca,'FontName','Arial','FontSize',12)
xlabel('Unchosen')
ylabel('Chosen')

%%% model-predicted confidence
%%% plot absdQ+Qc model
subplot(3,3,3)
idcModel = find(reg.whichLearnModel==ilearnmodel & ~reg.isConfPrev & strcmp(reg.confBias,'+Qc')); %sigmaQ
dum = squeeze(reg.confmodelmatChoice(idcModel,:,:,:))*100;
dum(dum==0) = NaN;
b = imagesc(squeeze(nanmean(dum))');
set(b,'AlphaData',~isnan(squeeze(nanmean(dum))));
set(b,'AlphaData',~eye(size(myMatMean)));

xticklabels(lg);yticklabels(lg);xtickangle(90)
xticks(1:8);yticks(1:8)
set(gca,'FontName','Arial','FontSize',12)
xlabel('Unchosen')
ylabel('Chosen')

caxis([50,90]);
h = colorbar()
% title(h,'\SigmaQ')
xticklabels(lg);yticklabels(lg);xtickangle(90)
xticks(1:8);yticks(1:8)
set(gca,'FontName','Arial','FontSize',12)
xlabel('Unchosen')
ylabel('Chosen')

subplot(3,3,6)
idcModel = find(reg.whichLearnModel==ilearnmodel & ~reg.isConfPrev & strcmp(reg.confBias,'+QcplusQu')); %sigmaQ
dum = squeeze(reg.confmodelmatChoice(idcModel,:,:,:))*100;
dum(dum==0) = NaN;
b = imagesc(squeeze(nanmean(dum))');
% set(b,'AlphaData',~isnan(dum));
set(b,'AlphaData',~eye(size(myMatMean)));

caxis([50,90]);
h = colorbar()
% title(h,'\SigmaQ')
xticklabels(lg);yticklabels(lg);xtickangle(90)
xticks(1:8);yticks(1:8)
set(gca,'FontName','Arial','FontSize',12)
xlabel('Unchosen')
ylabel('Chosen')

subplot(3,3,9)
idcModel = find(reg.whichLearnModel==ilearnmodel & ~reg.isConfPrev & strcmp(reg.confBias,'+V')); %sigmaQ
dum = squeeze(reg.confmodelmatChoice(idcModel,:,:,:))*100;
dum(dum==0) = NaN;
b = imagesc(squeeze(nanmean(dum))');
% set(b,'AlphaData',~isnan(dum));
set(b,'AlphaData',~eye(size(myMatMean)));

caxis([50,90]);
h = colorbar()
% title(h,'\SigmaQ')
xticklabels(lg);yticklabels(lg);xtickangle(90)
xticks(1:8);yticks(1:8)
set(gca,'FontName','Arial','FontSize',12)
xlabel('Unchosen')
ylabel('Chosen')

set(gcf,'Position',[438          68        1151         910])
% saveas(gcf,'Plots\hiddenVarsTransfer.svg')

%% differences between models

figure()
biasmodels = {7,9,10}; 
labels = {'Qc - unbiased','Qc - sigmaQ','Qc-V'};

dumqc = squeeze(reg.confmodelmatChoice(2,:,:,:))*100;

for i = 1:3
    subplot(3,1,i)
    dumi = squeeze(reg.confmodelmatChoice(biasmodels{i},:,:,:))*100;
    dum = dumqc-dumi;
    
    dum(dum==0) = NaN;
    b = imagesc(squeeze(nanmean(dum))');
%     set(b,'AlphaData',~isnan(dum));
    set(b,'AlphaData',~eye(8,8));
    
%     imagesc(squeeze(nanmean(dum))')
%     title(labels{i})
    
    dumavg = nanmean(dum);
    caxis([-1,1]*max(dumavg(:)));
    
    h = colorbar()
    title(h,labels{i})
    xticklabels(lg);yticklabels(lg);xtickangle(90)
    xticks(1:8);yticks(1:8)
    set(gca,'FontName','Arial','FontSize',12)
    xlabel('Unchosen')
    ylabel('Chosen')
end
    set(gcf,'Position', [897    42   343   954])
    
% saveas(gcf,'Plots\confmodelsTransferDifference.svg')

