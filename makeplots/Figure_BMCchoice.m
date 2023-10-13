%%% plot results of Bayesian Model Comparison for choice, for original
%%% model space
%%% 1. Comparison of all 
%%% 2. Family-wise comparison for REL vs no-REL, CONF vs no-CONF, and
%%% 3. ABS vs REL vs CONF vs RELCONF
%%% 4. variants within RELCONF

clear
subsel = 1:207;

learnTypes = {'ABS','CONF','REL','CONFREL'};
RuTypes = {'_0','_{Qu}','_{Rother}','_{Rlast}'};

%% load data

rlvars_all= load('Results\RLvars_all.mat');
LAME_all = [rlvars_all.LAME_allexp(:,[1:5,8:12])];
%% compare all learning models
options = struct();
options.DisplayWin = false;
[postBMC,outBMC] = VBA_groupBMC(-LAME_all(subsel,:)',options);
modelnames = {'ABS','REL_0','REL_{Qu}','REL_{ROther}','REL_{RLast}','CONF','CONFREL_0','CONFREL_{Qu}','CONFREL_{ROther}','CONFREL_{RLast}'};

figure()
set(gcf,'Color',[1,1,1])

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
    'XTickLabels',modelnames,...
    'XLim',[0 size(outBMC.Ef,1)+1],...
    'FontName','Arial')
xtickangle(90)
hL = legend('Exceedence P.','Protected EP','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff

saveas(gcf,'Plots/BMCchoice_Reduced.svg')
%% family-wise comparison REL vs no REL
options = struct();
options.families ={[1,6],[2:5,7:10]};
% options.families ={[1,8],[2:7,9:14]};

options.DisplayWin = false;
[postBMC,outBMC] = VBA_groupBMC(-LAME_all(subsel,:)',options);

famNames = {'no REL','REL'};
figure()
set(gcf,'Color',[1,1,1])
subplot(1,3,1)

hold on

bar(squeeze(100*outBMC.families.ep),...
    'FaceColor',.7*[1,1,1],...
    'EdgeColor','none')
alpha(0.5)

errorbar(100*outBMC.families.Ef,100*sqrt(outBMC.families.Vf),'d',...
    'Color',.5*[1,1,1],...
    'MarkerFaceColor',.7*[1,1,1],...
    'MarkerEdgeColor',.5*[1,1,1],...
    'LineStyle','None')
plot([-1,numel(outBMC.families.Ef)+1],100*[1./numel(outBMC.families.Ef),1./numel(outBMC.families.Ef)],'r--');
plot([-1,size(outBMC.families.Ef,1)+1],100*[0.95,0.95],'b--')

set(gca,'YLim', [0 100],...
    'XTick',1:numel(outBMC.families.Ef),...
    'XTickLabels',famNames,...
    'XLim',[0 size(outBMC.families.Ef,1)+1],...
    'FontName','Arial')
xtickangle(90)

hL = legend('Exceedence P.','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff

%% family-wise comparison CONF vs noCONF
options = struct();
% options.families ={[1:7],[8:14]};
options.families ={[1:5],[6:10]};

options.DisplayWin = false;
[postBMC,outBMC] = VBA_groupBMC(-LAME_all(subsel,:)',options);

famNames = {'no CONF','CONF'};
% figure()
subplot(1,3,2)
hold on

bar(squeeze(100*outBMC.families.ep),...
    'FaceColor',.7*[1,1,1],...
    'EdgeColor','none')
alpha(0.5)

errorbar(100*outBMC.families.Ef,100*sqrt(outBMC.families.Vf),'d',...
    'Color',.5*[1,1,1],...
    'MarkerFaceColor',.7*[1,1,1],...
    'MarkerEdgeColor',.5*[1,1,1],...
    'LineStyle','None')
plot([-1,numel(outBMC.families.Ef)+1],100*[1./numel(outBMC.families.Ef),1./numel(outBMC.families.Ef)],'r--');
plot([-1,size(outBMC.families.Ef,1)+1],100*[0.95,0.95],'b--')

set(gca,'YLim', [0 100],...
    'XTick',1:numel(outBMC.families.Ef),...
    'XTickLabels',famNames,...
    'XLim',[0 size(outBMC.families.Ef,1)+1],...
    'FontName','Arial')
xtickangle(90)

hL = legend('Exceedence P.','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff

%% family per counterfactual contextualization implementation
a = [2:5,7:10];
options = struct();
options.families ={[1,5],[2,6],[3,7],[4,8]};
options.DisplayWin = false;
[postBMC,outBMC] = VBA_groupBMC(-LAME_all(subsel,a)',options);

famNames = {'0','Q_U','Other','Last'};
% figure()
subplot(1,3,3)

hold on

bar(squeeze(100*outBMC.families.ep),...
    'FaceColor',.7*[1,1,1],...
    'EdgeColor','none')
alpha(0.5)

errorbar(100*outBMC.families.Ef,100*sqrt(outBMC.families.Vf),'d',...
    'Color',.5*[1,1,1],...
    'MarkerFaceColor',.7*[1,1,1],...
    'MarkerEdgeColor',.5*[1,1,1],...
    'LineStyle','None')
plot([-1,numel(outBMC.families.Ef)+1],100*[1./numel(outBMC.families.Ef),1./numel(outBMC.families.Ef)],'r--');
plot([-1,size(outBMC.families.Ef,1)+1],100*[0.95,0.95],'b--')

set(gca,'YLim', [0 100],...
    'XTick',1:numel(outBMC.families.Ef),...
    'XTickLabels',famNames,...
    'XLim',[0 size(outBMC.families.Ef,1)+1],...
    'FontName','Arial')
xtickangle(90)

hL = legend('Exceedence P.','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff

saveas(gcf,'Plots/BMCchoiceFamilies1_Reduced.svg')
%%%%%%%%%%%%

%% ABS vs REL vs CONF vs RELCONF
options = struct();
options.families ={1,[2:5],6,[7:10]};
options.DisplayWin = false;
[postBMC,outBMC] = VBA_groupBMC(-LAME_all(subsel,:)',options);
figure()
set(gcf,'Color',[1,1,1])

hold on

bar(squeeze(100*outBMC.families.ep),...
    'FaceColor',.7*[1,1,1],...
    'EdgeColor','none')
alpha(0.5)

errorbar(100*outBMC.families.Ef,100*sqrt(outBMC.families.Vf),'d',...
    'Color',.5*[1,1,1],...
    'MarkerFaceColor',.7*[1,1,1],...
    'MarkerEdgeColor',.5*[1,1,1],...
    'LineStyle','None')
plot([-1,numel(outBMC.families.Ef)+1],100*[1./numel(outBMC.families.Ef),1./numel(outBMC.families.Ef)],'r--');
plot([-1,size(outBMC.families.Ef,1)+1],100*[0.95,0.95],'b--')

set(gca,'YLim', [0 100],...
    'XTick',1:numel(outBMC.families.Ef),...
    'XTickLabels',learnTypes,...
    'XLim',[0 size(outBMC.families.Ef,1)+1],...
    'FontName','Arial')
xtickangle(90)

hL = legend('Exceedence P.','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff

saveas(gcf,'Plots/BMCchoiceRELvariants_Reduced.svg')

%% compare within winning family
options = struct();
options.DisplayWin = false;
[postBMC,outBMC] = VBA_groupBMC(-LAME_all(subsel,7:10)',options);
figure()
set(gcf,'Color',[1,1,1])
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
    'XTickLabels',{modelnames{7:10}},...
    'XLim',[0 size(outBMC.Ef,1)+1],...
    'FontName','Arial')

xtickangle(90)

hL = legend('Exceedence P.','Protected EP','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff
saveas(gcf,'Plots/BMCchoiceFamilies2_Reduced.svg')
