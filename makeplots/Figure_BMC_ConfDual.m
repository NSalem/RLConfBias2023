%%% model comparison for dual confidence models (single models predicting
%%% confidence in both tasks).
%%% Does family-wise comparison for same vs different parameter, for each
%%% parameter (intercept, difficulty coefficient, valence bias coefficient)

clear all
% close all force
clc

%% comparison for dual model (all within Qc model of bias)
load('Results\reg_conf_dual.mat')
figure()
set(gcf,'Color',[1,1,1])
%% diff vs same B0 within Qc model
a = find(whichConfModelLT == 1 & whichConfModelTT ==1 & whichLearnModel == 11)
options = struct();
options.DisplayWin = false;
options.families = {find(~isDiffB0(a)),find(isDiffB0(a))};
[post,out] =  VBA_groupBMC(-BICCONF(:,a)'/2,options);

subplot(1,3,1)

hold on

bar(squeeze(100*out.families.ep),...
    'FaceColor',.7*[1,1,1],...
    'EdgeColor','none')
alpha(0.5)

errorbar(100*out.families.Ef,100*sqrt(out.families.Vf),'d',...
    'Color',.5*[1,1,1],...
    'MarkerFaceColor',.7*[1,1,1],...
    'MarkerEdgeColor',.5*[1,1,1],...
    'LineStyle','None')
plot([-1,numel(out.families.Ef)+1],100*[1./numel(out.families.Ef),1./numel(out.families.Ef)],'r--');
plot([-1,size(out.families.Ef,1)+1],100*[0.95,0.95],'b--')

set(gca,'YLim', [0 100],...
    'XTick',1:numel(out.families.Ef),...
    'XTickLabels',{'same B_0','different B_0'},...
    'XLim',[0 size(out.families.Ef,1)+1],...
    'FontName','Arial')

xtickangle(90)

hL = legend('Exceedence P.','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff

%% within winner (diff B0), diff vs same dQ
a = find(whichConfModelLT == 1 & whichConfModelTT ==1 & whichLearnModel == 11 & isDiffB0)
options = struct();
options.DisplayWin = false;
options.families = {find(~isDiffDQ(a)),find(isDiffDQ(a))};
[post,out] =  VBA_groupBMC(-BICCONF(:,a)'/2,options);

subplot(1,3,2)
hold on

bar(squeeze(100*out.families.ep),...
    'FaceColor',.7*[1,1,1],...
    'EdgeColor','none')
alpha(0.5)

errorbar(100*out.families.Ef,100*sqrt(out.families.Vf),'d',...
    'Color',.5*[1,1,1],...
    'MarkerFaceColor',.7*[1,1,1],...
    'MarkerEdgeColor',.5*[1,1,1],...
    'LineStyle','None')
plot([-1,numel(out.families.Ef)+1],100*[1./numel(out.families.Ef),1./numel(out.families.Ef)],'r--');
plot([-1,size(out.families.Ef,1)+1],100*[0.95,0.95],'b--')

set(gca,'YLim', [0 100],...
    'XTick',1:numel(out.families.Ef),...
    'XTickLabels',{'same B_{|\DeltaQ|}','different B_{|\DeltaQ|}'},...
    'XLim',[0 size(out.families.Ef,1)+1],...
    'FontName','Arial')

xtickangle(90)

hL = legend('Exceedence P.','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff

%% within winner (diff B0, same dQ), diff vs same Qc
a = find(whichConfModelLT == 1 & whichConfModelTT ==1 & whichLearnModel == 11 & isDiffB0 & ~isDiffDQ)
options = struct();
options.DisplayWin = false;
options.families = {find(~isDiffParams(a)),find(isDiffParams(a))};
[post,out] =  VBA_groupBMC(-BICCONF(:,a)'/2,options);

subplot(1,3,3)
hold on
bar(squeeze(100*out.families.ep),...
    'FaceColor',.7*[1,1,1],...
    'EdgeColor','none')
alpha(0.5)

errorbar(100*out.families.Ef,100*sqrt(out.families.Vf),'d',...
    'Color',.5*[1,1,1],...
    'MarkerFaceColor',.7*[1,1,1],...
    'MarkerEdgeColor',.5*[1,1,1],...
    'LineStyle','None')
plot([-1,numel(out.families.Ef)+1],100*[1./numel(out.families.Ef),1./numel(out.families.Ef)],'r--');
plot([-1,size(out.families.Ef,1)+1],100*[0.95,0.95],'b--')

set(gca,'YLim', [0 100],...
    'XTick',1:numel(out.families.Ef),...
    'XTickLabels',{'same B_{Q_C}','different B_{Q_C}'},...
    'XLim',[0 size(out.families.Ef,1)+1],...
    'FontName','Arial')

xtickangle(90)

hL = legend('Exceedence P.','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff

set(gcf,'Position',[0,0,800,400])
saveas(gcf,'Plots/BMCDual.svg')

% All same vs all different RELOTHER-Qc
idcSame = find(whichLearnModel==11 & ~isDiffParams & whichConfModelLT==1 & ....
    whichConfModelTT==1 & ~isDiffDQ &~isDiffB0);

idcDiff = find(whichLearnModel==11 & isDiffParams & whichConfModelLT==1 & ....
    whichConfModelTT==1 & isDiffDQ &isDiffB0);

options = struct();
options.DisplayWin = false;
[postBMC,outBMC] = VBA_groupBMC(-BICCONF(:,[idcSame,idcDiff])'/2,options); 

figure()
hold on
bar(squeeze(100*outBMC.ep),...
    'FaceColor',.7*[1,1,1],...
    'EdgeColor','none')
alpha(0.5)

plot(squeeze(100*outBMC.pxp),'-o',...
    'Color',.7*[1,1,1],...
    'MarkerFaceColor',.7*[1,1,1],...
    'MarkerEdgeColor',[1,1,1])

errorbar(100*outBMC.Ef,100*sqrt(outBMC.Vf),'d',...
    'Color',.5*[1,1,1],...
    'MarkerFaceColor',.7*[1,1,1],...
    'MarkerEdgeColor',.5*[1,1,1],...
    'LineStyle','None')
plot([-1,numel(outBMC.Ef)+1],100*[1./numel(outBMC.Ef),1./numel(outBMC.Ef)],'r--');
plot([-1,size(outBMC.Ef,1)+1],100*[0.95,0.95],'b--')

set(gca,'YLim', [0 100],...
    'XTick',1:numel(outBMC.Ef),...
    'XTickLabels',{'All same','All different'},...
    'XLim',[0 size(outBMC.Ef,1)+1],...
    'FontName','Arial')

xtickangle(90)

hL = legend('Exceedence P.','Protected EP','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff
saveas(gcf,'Plots/BMCDual2.svg')

clear