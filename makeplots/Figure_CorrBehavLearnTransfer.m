%%% correlate overconfidence and valence bias between tasks, show
%%% scatterplots


clear all

addpath('helperfuncs')

% close all force


% load('Results\data_all.mat');
% addpath('helperfuncs');
% 
% whichexp = 1:10;
% conf_mat_sess = [];
% prevconf_mat_sess = [];
% corr_mat_sess = [];
% pref = [];
% accPost = [];
% confpost =[];
% prefMat_all= [];
% confMat_all = [];
% accMat_all = [];
% for iexp =whichexp
%     whichsess = size(data_all(iexp).cho,2);
%     subjects    = data_all(iexp).subjects;           % 3 sess
%     n_sub       = length(subjects);
%     n_sess      = numel(whichsess);%size(data_all(iexp).con,2);
%     n_trials    = size([data_all(iexp).con],3);
%     expN        = data_all(iexp).expN;
%     
%     thisExp =  MakeMatrix_sess_cond_trl(data_all(iexp),1);
%     
%     
%     conf_mat_sess = [conf_mat_sess;thisExp.conf_mat_sess.*100];
%     prevconf_mat_sess = [prevconf_mat_sess;thisExp.prevconf_mat_sess.*100];
%     corr_mat_sess = [corr_mat_sess;thisExp.corr_mat_sess];
%     pref = [pref;thisExp.pref];
%     confpost = [confpost;thisExp.confpostcon];
%     accPost = [accPost;thisExp.accPost];
%     
%     %         [prefMat,accMat,confMat] = makePostMats(data_all(iexp).ss,data_all(iexp).aa,data_all(iexp).conf_post);
%     %          prefMat_all = [prefMat_all;prefMat];
%     %          accMat_all = [accMat_all;accMat];
%     %          confMat_all = [confMat_all;confMat];
% end


load('Results/RLvars_all.mat')

[prefMat_all,accMat_all,confMat_all,confMatPair_all] = makePostMats(ss,aa,confpost);
%
accPost = 100.*squeeze(nanmean(accMat_all./4,3));
confPost = 100.*squeeze(nanmean(confMatPair_all,3));   
accLT = squeeze(nanmean(nanmean(correct,2),4).*100);
confLT = squeeze(nanmean(nanmean(conf,2),4)).*100;
overconfLT = confLT-accLT;
overconfLT = overconfLT(:,[1,3,2,4]); %RP,PP,RC,PC;
pair_mat = [6 5;2 1;8 7;4 3;];
newIdc = [6 5 2 1 8 7 4 3];
confTT =  confPost;
confTT = confTT(:,newIdc);
accTT = accPost(:,newIdc);
overconfTT = confTT-accTT;
overconfTT = [overconfTT(:,newIdc)];

% newIdc = pair_mat';
% XTT = squeeze(nanmean(nanmean(conf_mat_sess,2),4))-squeeze(nanmean(nanmean(corr_mat_sess,2),4).*100);


%% AVERAGE overconfidence across tasks
%    corr(nanmean(xLT,2),nanmean(xTT,2))

h1 = figure('Units', 'pixels', ...
    'Color',[1,1,1]);
h2 = figure('Units', 'pixels', ...
    'Color',[1,1,1]);
h3 = figure('Units', 'pixels', ...
    'Color',[1,1,1]);
h4 = figure('Units', 'pixels', ...
    'Color',[1,1,1]);

%repeat for avg overconfidence and avg valence bias
for imeas = 1:4
    if imeas ==1
        meanL = nanmean(overconfLT,2);
        meanT = nanmean(overconfTT,2);
        measname = 'Overcnfidence';
        figure(1)
    elseif imeas ==2
        meanL = nanmean(confLT(:,[1,2])-confLT(:,[3,4]),2);
        meanT = nanmean(confTT(:,[3,4,7,8])-confTT(:,[1,2,5,6]),2);
        measname = 'ValenceBias';
        figure(2)
    elseif imeas ==3
        meanL = mean(accLT(1:90,:),2);
        meanT = mean(accTT(1:90,:),2);
        measname = 'Accuracy';
        figure(3)
    elseif imeas ==4
        meanL = mean(confLT,2);
        meanT = mean(confTT,2);
        measname = 'Confidence'
        figure(4)
    end
    
    subplot('position',[0.82,0.15,0.12,0.6])
    h = histogram(meanT,15,'Orientation','horizontal','EdgeColor',[1,1,1],'FaceColor',[0.7,0.7,0.7]);
    'BinLimits',[min([meanL;meanT]),max([meanL;meanT])],...
    set(gca,'Fontname','Arial',...
        'FontSize',10)
    ylim([min([meanL;meanT]),max([meanL;meanT])])
    yticks({})
    box off
    subplot('Position',[0.15,0.15,0.6,0.6])
    hold on
    %     title(paramnames{ip})
    plot([min([meanL;meanT]),max([meanL;meanT])],[min([meanL;meanT]),max([meanL;meanT])],'Color',[.5,.5,.5],'LineStyle','-')
    plot(meanL,meanT,'o',...
        'MarkerFaceColor',[1,1,1],...
        'MarkerEdgeColor',[0,0,0]);
    ylim([min([meanL;meanT]),max([meanL;meanT])])
    xlim([min([meanL;meanT]),max([meanL;meanT])])
    
    dumlm = fitlm(meanL,meanT,'RobustOpts','on');
    lmAll{imeas} = dumlm;
    dumX = [min([meanL;meanT]):0.1:max([meanL;meanT])]';
    [y,yci] = dumlm.predict(dumX);
    plot(dumX,y,'Color',[0.5,0.0,0.5]);
%     plot(dumX,yci,'LineStyle','--','Color',[0.5,0.0,0.5]);
    fill([dumX;flipud(dumX)],[yci(:,2);flipud(yci(:,1))],[0.5,0.0,0.5],'LineStyle','None','FaceAlpha',0.2);

    hX = xlabel([measname, ' Learning Task']);
    hY = ylabel([measname, ' Transfer Task']);
    set([hX,hY],'Fontname','Arial',...
        'FontSize',10)
    set(gca,'Fontname','Arial',...
        'FontSize',10)
    subplot('Position',[0.15,0.82,0.6,0.12])
    histogram(meanL,15,'EdgeColor',[1,1,1],'FaceColor', [0.7,0.7,0.7],...
        'BinLimits',[min([meanL;meanT]),max([meanL;meanT])]);
    set(gca,'Fontname','Arial',...
        'FontSize',10)
    xlim([min([meanL;meanT]),max([meanL;meanT])])
    xticks({})
    box off
    set(gcf,'Position',[100,100,400,400]);
    saveas(gcf,['Plots\corrTasks',measname,'.svg'])
end


%accuracy confidence data only

meanL = mean(accLT(1:90,:),2);
meanT = mean(accTT(1:90,:),2);
measname = 'Accuracy';
figure()


subplot('position',[0.82,0.15,0.12,0.6])
h = histogram(meanT,15,'Orientation','horizontal','EdgeColor',[1,1,1],'FaceColor',[0.7,0.7,0.7]);
'BinLimits',[min([meanL;meanT]),max([meanL;meanT])],...
set(gca,'Fontname','Arial',...
    'FontSize',10)
ylim([min([meanL;meanT]),max([meanL;meanT])])
yticks({})
box off
subplot('Position',[0.15,0.15,0.6,0.6])
hold on
%     title(paramnames{ip})
plot([min([meanL;meanT]),max([meanL;meanT])],[min([meanL;meanT]),max([meanL;meanT])],'Color',[.5,.5,.5],'LineStyle','-')
plot(meanL,meanT,'o',...
    'MarkerFaceColor',[1,1,1],...
    'MarkerEdgeColor',[0,0,0]);
ylim([min([meanL;meanT]),max([meanL;meanT])])
xlim([min([meanL;meanT]),max([meanL;meanT])])

dumlm = fitlm(meanL,meanT,'RobustOpts','on');
lmAllAcc = dumlm;
dumX = [min([meanL;meanT]):0.1:max([meanL;meanT])]';
[y,yci] = dumlm.predict(dumX);
plot(dumX,y,'Color',[0.5,0.0,0.5]);
%     plot(dumX,yci,'LineStyle','--','Color',[0.5,0.0,0.5]);
fill([dumX;flipud(dumX)],[yci(:,2);flipud(yci(:,1))],[0.5,0.0,0.5],'LineStyle','None','FaceAlpha',0.2);

hX = xlabel([measname, ' Learning Task']);
hY = ylabel([measname, ' Transfer Task']);
set([hX,hY],'Fontname','Arial',...
    'FontSize',10)
set(gca,'Fontname','Arial',...
    'FontSize',10)
subplot('Position',[0.15,0.82,0.6,0.12])
histogram(meanL,15,'EdgeColor',[1,1,1],'FaceColor', [0.7,0.7,0.7],...
    'BinLimits',[min([meanL;meanT]),max([meanL;meanT])]);
set(gca,'Fontname','Arial',...
    'FontSize',10)
xlim([min([meanL;meanT]),max([meanL;meanT])])
xticks({})
box off
set(gcf,'Position',[100,100,400,400]);
saveas(gcf,['Plots\corrTasks_Accuracy_ConfData.svg'])


rmpath('helperfuncs')