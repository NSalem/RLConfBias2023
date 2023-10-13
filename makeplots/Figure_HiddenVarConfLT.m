% close all
% clc
%%%Analyze selected model fit for all experiments %%%
%%% code for context-confidence paper

clear; %close all force
addpath('ModelingFuncs\');
resultsdir = ['Results', filesep];
datadir = ['Data',filesep];
addpath('helperfuncs');

whichLearnModels =[11]; %RELCONFother

plotconf = 1;
whichexp = [1:10];
subsel = [1:90];

% exps = load('Results\data_all')
% exps_conf = exps.data_all;
load('Results\RLvars_all.mat');
% load('Results\RLvars_norandom_all.mat');

% load('Results\RLvars_ExtraContext_all.mat');

reg = load('Results\reg_conflogit_learning_dqabs.mat');
% reg = load('Results\reg_conflogit_norandom_learning_dqabs.mat');

col_mat = [0,1,0;...
    0,1,1;...
    1,0,0;...
    1,0,1];
MS = 4;

Xfill = [1:size(correct,4),fliplr(1:size(correct,4))];


conf_mat_sess = reg.confmat;
prev_conf_mat = confprev;
corr_mat_sess = correct;
conf_mat    = squeeze(nanmean(reg.confmat(subsel,:,:,:),2));
conf_m = squeeze(nanmean(conf_mat));
conf_se = squeeze(nanstd(conf_mat)./sqrt(size(conf_mat,1)));


for nlearnmodel = whichLearnModels
    
    idc = [find(reg.whichLearnModel ==nlearnmodel  & reg.isConfPrev & ~(strcmp(reg.confBias,'')))]; %REL_Qc model, REL_sigmaQ,REL_V respectively
    confcoeff = reg.confcoeff;
    confmodels =idc;
    
    %     for iconfmodel = 1:numel(idc)
    %         nconfmodel = confmodels(iconfmodel);
    %         for isub = 1:numel(subsel)
    %             coeffsSub = confcoeff{nconfmodel}(isub,:);
    %             dQSub = squeeze(abs(dQ(nlearnmodel,isub,:,:,:)));
    %             if iconfmodel == 1
    %                 contextSub = squeeze(Q_c(nlearnmodel,isub,:,:,:));
    %             elseif iconfmodel == 2
    %                 contextSub = squeeze(Q_c(nlearnmodel,isub,:,:,:)+Q_uc(nlearnmodel,isub,:,:,:));
    %             else
    %                 contextSub = squeeze(V(nlearnmodel,isub,:,:,:));
    %             end
    %             prevconfSub = squeeze(prevconf_mat_sess(isub,:,:,:));
    %             y(iconfmodel,isub,:,:,:) = coeffsSub(1)+coeffsSub(2).*dQSub+coeffsSub(3).*contextSub+coeffsSub(4).*prevconfSub;
    %         end
    %     end
    
    corr_mat    = squeeze(nanmean(corr_mat_sess(subsel,:,:,:),2));
    pc_mat      = squeeze(nanmean(pc(nlearnmodel,subsel,:,:,:),3));
    
    %% plot hidden variables (Qc,Quc,sigma Q,dQ, and V)
    Qc_all = squeeze(nanmean(Q_c(nlearnmodel,:,:,:,:),3));
    Qc_m = squeeze(nanmean(Qc_all));
    Qc_se = squeeze(nanstd(Qc_all)./sqrt(size(Qc_all,1)));
    
    Qu_all = squeeze(nanmean(Q_uc(nlearnmodel,:,:,:,:),3));;
    Qu_m = squeeze(nanmean(Qu_all));
    Qu_se = squeeze(nanstd(Qu_all)./sqrt(size(Qu_all,1)));
    
    sigmaQ_all = Qc_all+Qu_all;
    sigmaQ_m = squeeze(nanmean(sigmaQ_all));
    sigmaQ_se = squeeze(nanstd(sigmaQ_all)./sqrt(size(sigmaQ_all,1)));
    
    %%% Qc minus Qu %%%
    QCMU = Q_c-Q_uc;
    QCMU_all = squeeze(nanmean(QCMU(nlearnmodel,:,:,:,:),3));
    QCMU_m = squeeze(nanmean(QCMU_all));
    QCMU_se = squeeze(nanstd(QCMU_all)./sqrt(size(QCMU_all,1)));
    %%%abs dQ%%%
    dQAbs_all = squeeze(nanmean(abs(dQ(nlearnmodel,:,:,:,:)),3));
    dQAbs_m = squeeze(nanmean(dQAbs_all));
    dQAbs_se = squeeze(nanstd(dQAbs_all)./sqrt(size(dQAbs_all,1)));
    
    %%% Best minus worst
    QBMW = squeeze(Q(:,:,:,:,2,:)  - Q(:,:,:,:,1,:));
    QBMW_all =squeeze(nanmean(QBMW(nlearnmodel,:,:,:,:),3));
    QBMW_m = squeeze(nanmean(QBMW_all));
    QBMW_se = squeeze(nanstd(QBMW_all)./sqrt(size(QBMW_all,1)));
    
    QMaxMin = squeeze(max(Q,[],5)-min(Q,[],5));
    QMaxMin_all =squeeze(nanmean(QMaxMin(nlearnmodel,:,:,:,:),3));
    QMaxMin_m = squeeze(nanmean(QMaxMin_all));
    QMaxMin_se = squeeze(nanstd(QMaxMin_all)./sqrt(size(QMaxMin_all,1)));
    
    V_all = squeeze(nanmean(V(nlearnmodel,:,:,:,:),3));
    V_m = squeeze(nanmean(V_all));
    V_se = squeeze(nanstd(V_all)./sqrt(size(V_all,1)));
    
    pc_m  = squeeze(nanmean(pc_mat));
    pc_se = squeeze(nanstd(pc_mat)./sqrt(size(pc_mat,1)));
    
    %% plot Qc, sigmaQ, V
    h1 = figure('Units', 'pixels', ...
        'Position', [400 200 1300 800],...
        'Color',[1,1,1]);
    %% dQ
    subplot(3,3,1)
    for k = [1:4]
        hold on
        fill(Xfill,[dQAbs_m(k,:)-dQAbs_se(k,:),fliplr(dQAbs_m(k,:)+dQAbs_se(k,:))],.65*col_mat(k,:),'EdgeColor','none')
        alpha(0.5)
    end
    hX = xlabel('Trials');
    hY = ylabel('|\Delta Q|');
    set(gca,'FontName','Arial','FontSize',10,'XLim',[1,24],'YLim', [0,0.75])
    set([hX,hY],'FontName','Arial','FontSize',10)
    subplot(3,3,2)
    for k = [1:4]
        hold on
        fill(Xfill,[Qc_m(k,:)-Qc_se(k,:),fliplr(Qc_m(k,:)+Qc_se(k,:))],.65*col_mat(k,:),'EdgeColor','none')
        alpha(0.5)
    end
    set(gca,'FontName','Arial','FontSize',10,'XLim',[1,24], 'YLim', [-0.25,0.5])
    hX = xlabel('Trials');
    hY = ylabel('Q_{c}');
    set([hX,hY],'FontName','Arial','FontSize',10)
    
    %sigmaQ
    subplot(3,3,5)
    for k = [1:4]
        hold on
        fill(Xfill,[sigmaQ_m(k,:)-sigmaQ_se(k,:),fliplr(sigmaQ_m(k,:)+sigmaQ_se(k,:))],.65*col_mat(k,:),'EdgeColor','none')
        alpha(0.5)
        
    end
    set(gca,'FontName','Arial','FontSize',10,'XLim',[1,24],'YLim', [-0.5,0.75])
    hX = xlabel('Trials');
    hY = ylabel('SigmaQ');
    set([hX,hY],'FontName','Arial','FontSize',10)
    
    %V
    subplot(3,3,8)
    for k = [1:4]
        hold on
        fill(Xfill,[V_m(k,:)-V_se(k,:),fliplr(V_m(k,:)+V_se(k,:))],.65*col_mat(k,:),'EdgeColor','none')
        alpha(0.5)
    end
    
    set(gca,'FontName','Arial','FontSize',10,'XLim',[1,24], 'YLim', [-.5,.5])
    hX = xlabel('Trials');
    hY = ylabel('V');
    set([hX,hY],'FontName','Arial','FontSize',10)
    
    
    %% confidence models
    if plotconf
        y = reg.confmodel(:,idc,:,:,:)*100;
        y = permute(y,[2,1,3,4,5]);
        %% V model
        if numel(idc)>2
            subplot(3,3,9)
            if numel(idc)>2
                dum  = squeeze((nanmean(y(3,:,:,:,:),3)));
                dum_m = squeeze(nanmean(dum));
                dum_se = squeeze(nanstd(dum)./sqrt(size(dum,1)));
                
                for k = [1:4]
                    hold on
                    fill(Xfill,[dum_m(k,:)-dum_se(k,:),fliplr(dum_m(k,:)+dum_se(k,:))],.65*col_mat(k,:),'EdgeColor','none')
                    alpha(0.5)
                    errorbar(conf_m(k,:),conf_se(k,:),'.','Color',.65*col_mat(k,:),...
                        'LineStyle','none');
                end
                set(gca,'FontName','Arial','FontSize',10,'XLim',[1,24], 'YLim', [50,100],'YTick',[50 75 100])
                hX = xlabel('Trials');
                set([hX],'FontName','Arial','FontSize',10)
            else
                hold on
                plot([0,24],[50,100],'k')
                plot([0,24],[100,50],'k')
                yticks([]); xticks([])
            end
        end
        
        dum  = squeeze((nanmean(y(1,:,:,:,:),3)));
        dum_m = squeeze(nanmean(dum));
        dum_se = squeeze(nanstd(dum)./sqrt(size(dum,1)));
        
        subplot(3,3,3)
        for k = [1:4]
            hold on
            fill(Xfill,[dum_m(k,:)-dum_se(k,:),fliplr(dum_m(k,:)+dum_se(k,:))],.65*col_mat(k,:),'EdgeColor','none')
            alpha(0.5)
            errorbar(conf_m(k,:),conf_se(k,:),'.','Color',.65*col_mat(k,:),...
                'LineStyle','none');
        end
        set(gca,'FontName','Arial','FontSize',10,'XLim',[1,24], 'YLim', [50,100],'YTick',[50 75 100])
        hX = xlabel('Trials');
        set([hX],'FontName','Arial','FontSize',10)
        
        dum  = squeeze((nanmean(y(2,:,:,:,:),3)));
        dum_m = squeeze(nanmean(dum));
        dum_se = squeeze(nanstd(dum)./sqrt(size(dum,1)));
        
        subplot(3,3,6)
        for k = [1:4]
            hold on
            fill(Xfill,[dum_m(k,:)-dum_se(k,:),fliplr(dum_m(k,:)+dum_se(k,:))],.65*col_mat(k,:),'EdgeColor','none')
            alpha(0.5)
            errorbar(conf_m(k,:),conf_se(k,:),'.','Color',.65*col_mat(k,:),...
                'LineStyle','none');
        end
        set(gca,'FontName','Arial','FontSize',10,'XLim',[1,24], 'YLim', [50,100],'YTick',[50 75 100])
        hX = xlabel('Trials');
        set([hX],'FontName','Arial','FontSize',10)
        set(gcf,'Position',[0,0 ,600,600])
    end
end

saveas(gcf,'plots/hiddenVarsLearn.svg')

%% plot deltaQ and pc side by side

figure()
for ip = 1:4
    subplot(1,4,ip)
    switch ip
        case 1
            dum_m = dQAbs_m;
            dum_se= dQAbs_se;
            yl = [0,1];
            lab = '|\DeltaQ|';
        case 2
            dum_m = QBMW_m;
            dum_se= QBMW_se;
            yl = [0,1];
            lab = 'Q_{BEST}-Q_{WORST}';
            
        case 3
            dum_m = QCMU_m;
            dum_se= QCMU_se;
            yl = [0,1];
            lab = 'Q_{C}-Q_{U}';
        case 4
            
            dum_m = pc_m;
            dum_se= pc_se;
            yl = [0.5,1];
            lab = 'p(correct)';
    end
    for k = [1:4]
        hold on
        fill(Xfill,[dum_m(k,:)-dum_se(k,:),fliplr(dum_m(k,:)+dum_se(k,:))],.65*col_mat(k,:),'EdgeColor','none')
        alpha(0.5)
        plot(dum_m(k,:), 'Color',col_mat(k,:))
        %         errorbar(conf_m(k,:),conf_se(k,:),'.','Color',.65*col_mat(k,:),...
        %             'LineStyle','none');
    end
    set(gca,'FontName','Arial','FontSize',12,'XLim',[1,24], 'YLim', yl,'YTick',[50 75 100])
    hX = xlabel('Trials');
    ylabel(lab)
    
end
saveas(gcf,'Plots/hiddenVarsDifficultyLearn.svg')

% set(gcf, 'Position',[0,0,1500,400])
%
% saveas(gcf,'Plots/deltaQvsPC.svg');

%% regressions

if plotconf
    %regress avg across sesions
    myX = squeeze(nanmean(nanmean(squeeze(y(1,:,:,:,:)),1),2));
    myY = squeeze(nanmean(nanmean(conf_mat_sess,1),2)).*100;
    dumlm = fitlm(myX(:),myY(:),'RobustOpts','on')


    %regress w/o averaging across sessions
    myX = squeeze(nanmean(squeeze(y(1,:,:,:,:)),1));
    myY = squeeze(nanmean(conf_mat_sess,1)).*100;
    dumlm2 = fitlm(myX(:),myY(:),'RobustOpts','on')


    for isub = 1:size(y,2)
        myY = squeeze(nanmean(y(1,isub,:,:,:),3))*100;
        myX = squeeze(nanmean(conf_mat_sess(isub,:,:,:),2));
        rall(isub) = corr(myX(:),myY(:));
    end
    
end


rmpath('helperfuncs');

