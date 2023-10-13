%%% plot conf parameter distributions for Qc+Qu models %%%

addpath('helperfuncs');

%% confidence model params
regLT = load('Results\reg_conflogit_learning_nodq.mat'); %learning task
idc =10;
coeffLT = regLT.confcoeff{idc};
% coeffLT(:,5) = coeffLT(:,5)*100;

regTT = load('Results\reg_conflogit_posttest_nodq.mat'); %transfer task
idc = 5;
coeffTT = regTT.confcoeff{idc};
% coeffTT(:,4) = coeffTT(:,4)*100;

col = repmat(.5,5,3);
figure()
subplot(1,2,1)
hold on;
plot([0,5],[0,0],'k--')
pirateplot(coeffLT',col,-15,110,12,'Learning Task','','');
xticklabels({'B_{0}','B_{Q_C}','B_{Q_U}','B_V','B_{Conf(t-1)}'});
xtickangle(90);
ylim([-20,20])
subplot(1,2,2)
hold on
plot([0,5],[0,0],'k--')
pirateplot(coeffTT',col,-15,110,12,'Transfer Task','','');
xticklabels({'B_{0}','B_{Q_C}','B_{Q_U}','B_V','B_{Conf(t-1)}'});
xtickangle(90);
ylim([-20,20])

saveas(gcf,'Plots/paramsconfQcQuV.svg');

%%% t-tests 
% all params vs 0
[~,pL,ciL,tL] =  ttest(coeffLT);
[~,pT,ciT,tT] =  ttest(coeffTT);

% Qc vs -Qu
[~,pCU_L,ciCU_L,tCU_L] =  ttest(coeffLT(:,2),-coeffLT(:,3));
[~,pCU_T,ciCU_T,tCU_T] =   ttest(coeffTT(:,2),-coeffTT(:,3));

%% Qc+Qu

%%% confidence model params
regLT = load('Results\reg_conflogit_learning_nodq.mat'); %learning task
idc =9;
coeffLT = regLT.confcoeff{idc};
% coeffLT(:,4) = coeffLT(:,4)*100;

regTT = load('Results\reg_conflogit_posttest_nodq.mat'); %transfer task
idc = 4;
coeffTT = regTT.confcoeff{idc};
% coeffTT(:,4) = coeffTT(:,4)*100;

col = repmat(.5,5,3);
figure()
subplot(1,2,1)
hold on;
plot([0,5],[0,0],'k--')
pirateplot(coeffLT',col,-15,110,12,'Learning Task','','');
xticklabels({'B_{0}','B_{Q_C}','B_{Q_U}','B_{Conf(t-1)}'});
xtickangle(90);
ylim([-20,20])
subplot(1,2,2)
hold on
plot([0,5],[0,0],'k--')
pirateplot(coeffTT',col,-15,110,12,'Transfer Task','','');
xticklabels({'B_{0}','B_{Q_C}','B_{Q_U}','B_{Conf(t-1)}'});
xtickangle(90);
ylim([-20,20])

saveas(gcf,'Plots/paramsconfQcQu.svg');

% rmpath('helperfuncs')
