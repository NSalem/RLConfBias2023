function plotModelFitsLearn(rlvarsfile, confregfile, subsel)

if nargin <2 || isempty(confregfile)
    confregfile = 'Results/reg_conflogit_learning_dqabs.mat';
end

if nargin<1 || isempty (rlvarsfile)
    rlvarsfile = 'Results/RLvars_all.mat';
end
%% Load data
% load('Results\data_all.mat');
addpath('helperfuncs');

rlVars = load(rlvarsfile);
reg = load(confregfile);

if nargin<3 || isempty(subsel)
    subsel = 1:size(rlVars.correct,1);
end

imodel = 11;
pcMod = squeeze(nanmean(rlVars.pc(imodel,subsel,:,:,:),3)).*100;
% clear rlVars;

idc = find(reg.isConfPrev & reg.whichLearnModel==imodel & strcmp(reg.confBias, '+Qc'));
confMod = squeeze(nanmean(reg.confmodel(subsel,idc,:,:,:),3));
clear reg
% whichexp = 1:10;
% conf_mat_sess = [];
% corr_mat_sess = [];
% for iexp =whichexp
%     whichsess = 1:size(data_all(iexp).corr,2);
%     subjects    = data_all(iexp).subjects;           % 3 sess
%     n_sub       = length(subjects);
%     n_sess      = numel(whichsess);%size(data_all(iexp).con,2);
%     n_trials    = size([data_all(iexp).con],3);
%     expN        = data_all(iexp).expN;
%     
%     thisExp =  MakeMatrix_sess_cond_trl(data_all(iexp),1);
%     conf_mat_sess = [conf_mat_sess;thisExp.conf_mat_sess.*100];
%     corr_mat_sess = [corr_mat_sess;thisExp.corr_mat_sess];
% end

%% Average session data
% corr_mat    = squeeze(nanmean(corr_mat_sess,2));
% conf_mat    = squeeze(nanmean(conf_mat_sess,2));

corr_mat = squeeze(nanmean(rlVars.correct(subsel,:,:,:),2));
conf_mat =  squeeze(nanmean(rlVars.conf(subsel,:,:,:),2));

%% Plot results
h1 = figure('Units', 'pixels', ...
    'Position', [300 100 800 700],...
    'Color',[1,1,1]);

col_mat = .65*[0,1,0;...
    0,1,1;...
    1,0,0;...
    1,0,1];

reorderX = [1 3 2 4];
ebarpos = .3*[-1 +1 -1 +1];
Xfill = [1:24,fliplr(1:24)];

for k_meas = [1:3]
    switch k_meas
        case 1
            X       = 100*corr_mat(:,:,:);
            Xmod    = pcMod(:,:,:);
            yTTL    = 'Accuracy (%)';
            yl = [40 100];
            yl2 = [20 100];
            %             plot([1 24],[.5,.5],':k')
        case 2
            X       = 100*conf_mat(:,:,:);
            Xmod    = confMod*100;
            yTTL = 'Confidence (%)';
            yl = [40 100];
            yl2 = [50 100];
        case 3
            X       = 100*conf_mat(:,:,:)-100*corr_mat(:,:,:);
            Xmod    = 100*confMod-pcMod(:,:,:);
            yTTL = 'Overconfidence (%)';
            yl = [-20 30];
            yl2 = [-50 75];
    end
    
    Xt_m      = squeeze(nanmean(X));
    Xt_se     = squeeze(nanstd(X))./sqrt(size(X,1));


    %%% stats %%%
    Xx = squeeze(nanmean(X(:,:,1:24),3)); %avg across time
    XxMod = squeeze(nanmean(Xmod(:,:,1:24),3));

    n_sub = size(Xx,1);
    vnames = {'GainLoss','Completepartial','subject'};
    f1 = [zeros(2*n_sub,1);ones(2*n_sub,1)];
    f2 = repmat([zeros(n_sub,1);ones(n_sub,1)],2,1);
    rfx = (1:n_sub)';
    f3 = repmat(rfx,4,1);
    [pdata,anovatabledata{k_meas},stats]=anovan(Xx(:),{f1,f2,f3},'random',3,'model','interaction','varnames',vnames,'display', 'on');
    [cdata,m,h,nms] = multcompare(stats,'Dimension',[1,2],'display', 'off');
    congroupletters{k_meas} = makeGroupsMultipleComparisons(cdata);
    Fidx = [2,3,5];
    for iF = [1,2,3]
        nF = Fidx(iF);
        F = anovatabledata{k_meas}{nF,6};
        pF =  anovatabledata{k_meas}{nF,7};
        myF(k_meas,iF) = F;
        myP(k_meas,iF) = pF;
    end
    
    
    %%% plot partial conditions
    subplot(3,3,sub2ind([3,3],1,k_meas))
    hold on
    for k = [1 3]
        errorbar(Xt_m(k,:),Xt_se(k,:),'.','Color',.65*col_mat(k,:),...
                    'LineStyle','none');
        mtp = squeeze(nanmean(Xmod(:,k,:),1));
        stp = squeeze(nanstd(Xmod(:,k,:),0,1))./sqrt(size(Xmod,1));
        fill(Xfill,[mtp'-stp',fliplr(mtp'+stp')],.65*col_mat(k,:),'EdgeColor','none')
        alpha(0.5)
%         plot(mtp,'Color',col_mat(k,:));
    end
    set(gca,'XLim',[0 25],...
        'YLim',yl,...
        'YTick',min(yl):10:max(yl),...
        'FontName','Arial',...
        'FontSize',10)
    title('Partial')
    hX = xlabel('Trials');
    hY = ylabel(yTTL);
    set([hX,hY],...
        'FontName','Arial',...
        'FontSize',11)
    
    if k_meas ==3
        plot([0,25],[0,0],'k--')
    else
       plot([0,25],[50,50],'k--') 
    end
    
    %%% plot complete conditions
    subplot(3,3,sub2ind([3,3],2,k_meas))
    hold on
    for k = [2 4];
        errorbar(Xt_m(k,:),Xt_se(k,:),'.','Color',.65*col_mat(k,:),...
                'LineStyle','none');
        mtp = squeeze(nanmean(Xmod(:,k,:),1));
        stp = squeeze(nanstd(Xmod(:,k,:),0,1))./sqrt(size(Xmod,1));
        fill(Xfill,[mtp'-stp',fliplr(mtp'+stp')],.65*col_mat(k,:),'EdgeColor','none')
        alpha(0.5)
%         plot(mtp,'Color',col_mat(k,:));
    end
    set(gca,'XLim',[0 25],...
        'YLim',yl,...
        'YTick',min(yl):10:max(yl),...
        'FontName','Arial',...
        'FontSize',10)
    title('Complete')
    hX = xlabel('Trials');
    hY = ylabel(yTTL);
    set([hX,hY],...
        'FontName','Arial',...
        'FontSize',11)
    
    if k_meas ==3
        plot([0,25],[0,0],'k--')
        else
       plot([0,25],[50,50],'k--') 
    end
    
    %%% plot violins
    subplot(3,3,sub2ind([3,3],3,k_meas));
    dum = squeeze(nanmean(Xx,3));
    dotproperties.useJitter = 1;
    dotproperties.plotLines = [1,2;3,4];
    pirateplot(dum(:,[1,3,2,4])',col_mat([1,3,2,4],:),yl(1),yl(2),12,'','', yTTL,dotproperties);
    errorbar([1,3,2,4],nanmean(XxMod),std(XxMod)./sqrt(size(XxMod,1)),...
        'LineStyle','none','Marker','o','MarkerFaceColor','white',...
        'MarkerEdgeColor','black','Color','black');
    
    if k_meas ==3
        plot([0,5],[0,0],'k--')
       else
       plot([0,5],[50,50],'k--') 
    end
    %     text(0.1,-.05,spr1intf('p_{Valence} = %0.3f',myP(k_meas,1)),'Units', 'normalized');
    %     text(0.1,-.15,sprintf('p_{Information} = %0.3f',myP(k_meas,2)),'Units', 'normalized');
    %     text(0.1,-.25,sprintf('p_{Interaction} = %0.3f',myP(k_meas,3)),'Units', 'normalized');
%     
%     for icond = 1:4
%         text(icond,yl(2)-.05*(yl(2)-yl(1)),congroupletters{k_meas}{icond},'HorizontalAlignment','center');
%     end
    
    meanAll = nanmean(nanmean(dum,2),1);
    seAll = nanstd(nanmean(dum,2),1)./sqrt(size(dum,1));
    errorbar(5,meanAll,seAll,'Color','k','Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k')
    xlim([0,6])
    
    set(gca,...
        'Ylim',yl2,...
        'FontName','Arial',...
        'FontSize',10)
    
end


% saveas(gcf,'Plots/modelfitsLearn.svg');


end