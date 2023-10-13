%%% correlations and differences in confidence parameters for learning and
%%% transfer task (models fitted separatelty for each task.
%%% note that the model comparison in Figure10_BMCconf_TTLT.m is for the
%%% "dual" models that are fitted jointly on both tasks)

clear
% close all force
clc

addpath('helperfuncs\')

confLearn = load('Results\reg_conflogit_learning_dqabs.mat');
confPost = load('Results\reg_conflogit_posttest_dqabs.mat');

learnomdoel = 11;
confmodelL = "conf ~ 1 + dQabs+Qc+ confprev";
confmodelidcL = strcmp(confLearn.formulas,confmodelL);
bLearn = confLearn.confcoeff{confmodelidcL};
tLearn = confLearn.conft{confmodelidcL};
bLearn = bLearn(:,1:3); 

confmodelT = "conf ~ 1 + dQabs+Qc";
confmodelidcT = strcmp(confPost.formulas,confmodelT);
bPost = confPost.confcoeff{confmodelidcT};
tPost = confPost.conft{confmodelidcT};

paramnames = {'b_0','b_{|\DeltaQ|}','b_{Q_C}'};

dotproperties.useJitter = 1;

for ip = 1:size(bLearn,2)
    figure()
    subplot('position',[0.82,0.15,0.12,0.6])
    h = histogram(bPost(:,ip),15,'Orientation','horizontal','EdgeColor',[1,1,1],'FaceColor', [0.7,0.7,0.7],...
        'BinLimits',[min([bLearn(:,ip);bPost(:,ip)]),max([bLearn(:,ip);bPost(:,ip)])]);
    ylim([min([bLearn(:,ip);bPost(:,ip)]),max([bLearn(:,ip);bPost(:,ip)])])
    yticks({})
    box off
    subplot('Position',[0.15,0.15,0.6,0.6])
    hold on
    title(paramnames{ip})
    scatter(bLearn(:,ip),bPost(:,ip),'o',...
        'MarkerFaceColor',[1,1,1],...
        'MarkerEdgeColor',[0,0,0]);
    ylim([min([bLearn(:,ip);bPost(:,ip)]),max([bLearn(:,ip);bPost(:,ip)])])
    xlim([min([bLearn(:,ip);bPost(:,ip)]),max([bLearn(:,ip);bPost(:,ip)])])
    plot([min([bLearn(:,ip);bPost(:,ip)]),max([bLearn(:,ip);bPost(:,ip)])],[min([bLearn(:,ip);bPost(:,ip)]),max([bLearn(:,ip);bPost(:,ip)])],'Color',[.5,.5,.5],'LineStyle',':')
    dumlm = fitlm(bLearn(:,ip),bPost(:,ip),'RobustOpts','on');
    lmAll{ip} = dumlm;
    
    [corrAll(ip),corrPAll(ip)] = corr(bLearn(:,ip),bPost(:,ip));
    [~,tPAll(ip),~,tAll(ip)]=ttest(bPost(:,ip),bLearn(:,ip));
    dumX = [min([bLearn(:,ip);bPost(:,ip)]):0.1:max([bLearn(:,ip);bPost(:,ip)])]';
    [y,yci] = dumlm.predict(dumX);
    plot(dumX,y,'Color',[0.5,0.0,0.5]);
    %     plot(dumX,yci,'LineStyle','--','Color',[0.5,0.0,0.5]);
    patch([dumX;flipud(dumX)],[yci(:,2);flipud(yci(:,1))],[0.5,0.0,0.5],'LineStyle','None','FaceAlpha',0.2);
    xlabel('Learning task')
    ylabel('Transfer task')
    set(gca,'FontName','Arial',...
        'FontSize',16)
    subplot('Position',[0.15,0.82,0.6,0.12])
    histogram(bLearn(:,ip),15,'EdgeColor',[1,1,1],'FaceColor', [0.7,0.7,0.7],...
        'BinLimits',[min([bLearn(:,ip);bPost(:,ip)]),max([bLearn(:,ip);bPost(:,ip)])]);
    xlim([min([bLearn(:,ip);bPost(:,ip)]),max([bLearn(:,ip);bPost(:,ip)])])
    set(gca,'FontName','Arial')
    xticks({})
    box off
    set(gcf,'Position',[100,100,400,400],...
        'Color',[1,1,1]);
    saveas(gcf,['Plots/corrConfParamsLT_TT',num2str(ip),'.svg'])
end

% 
% figure()
% diffmax = max([bPost(:)-bLearn(:)])
% diffmin = min([bPost(:)-bLearn(:)])
% pirateplot(bPost'-bLearn',repmat([.5,.5,.5],3,1),diffmin,diffmax,14,'Transfer - Learning','','',dotproperties);
% xticklabels(paramnames);
% hold on
% plot([0,4],[0,0],'k--')
% set(gcf,'Position',[100,100,500,500]);

rmpath('helperfuncs')
