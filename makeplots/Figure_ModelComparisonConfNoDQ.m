

%% Learning task 

%% Learning task 
load('Results\reg_conflogit_learning_nodq.mat')
LL_L = LLCONF;
npars_L = npars;
BIC_L = BICCONF;

confmodelnames ={'Qc','\SigmaQ','V','Qc+Qu','Qc+Qu+V'};
%% quantitative model comparison 
options = struct();
options.DisplayWin = false;
idc = 6:10;
[postBMC,outBMC] = VBA_groupBMC(-BICCONF(:,idc)'/2,options);


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
    'XTickLabels',confmodelnames,...
    'XLim',[0 size(outBMC.Ef,1)+1],...
    'FontName','Arial')

xtickangle(90)

hL = legend('Exceedence P.','Protected EP','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff
title('learning')

%% Model comparison for transfer task 

clear LLCONF BICCONF npars
load('Results\reg_conflogit_posttest_nodq.mat')
LL_T = LLCONF;
BIC_T = BICCONF;
npars_T = npars;

confmodelnames ={'Qc','\SigmaQ','V','Qc+Qu','Qc+Qu+V'};
%% quantitative model comparison 
options = struct();
options.DisplayWin = false;
idc = 1:5;
[postBMC,outBMC] = VBA_groupBMC(-BICCONF(:,idc)'/2,options);

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
    'XTickLabels',confmodelnames,...
    'XLim',[0 size(outBMC.Ef,1)+1],...
    'FontName','Arial')

xtickangle(90)

hL = legend('Exceedence P.','Protected EP','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff
title('transfer')


%% joint BIC for Learning and Transfer
LLallJoint = LL_L+LL_T;
BICallJoint = (npars_L+npars_T)*log(112+288)-2*LLallJoint;

idc = 6:10; %use models with autocorrelation term
[postBMC,outBMC] = VBA_groupBMC(-BICallJoint(:,idc)'/2,options);

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
    'XTickLabels',confmodelnames,...
    'XLim',[0 size(outBMC.Ef,1)+1],...
    'FontName','Arial')

xtickangle(90)

hL = legend('Exceedence P.','Protected EP','Expected F.');
hY = ylabel('probability');
set([hY,hL],'FontName','Arial')
set(hL,'Location','Best')
legend boxoff
title('transfer')
