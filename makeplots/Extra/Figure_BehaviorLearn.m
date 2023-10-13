%Plot participant behavior in learning task across conditions (curves and 
%violin plots for accuracy and confidence)

% clc;
clear;
close all force;

addpath('helperfuncs')

%% Load data

load('Results\RLvars_all.mat');
exp_mat_sess = infostr.subexp'.*ones(size(correct));

%% Average across sessions 
corr_mat    = squeeze(nanmean(correct,2));
conf_mat    = squeeze(nanmean(conf,2));

%% Plot results
h1 = figure('Units', 'pixels', ...
    'Position', [400 200 800 700],...
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
            X       = 100*corr_mat;
            yTTL    = 'Accuracy (%)';
            yl = [40 100];
            yl2 = [20 100];
            f4 = categorical(squeeze(exp_mat_sess(:,1,:,1)));

            %             plot([1 24],[.5,.5],':k')
        case 2
            X       = conf_mat(1:90,:,:)*100;
            yTTL = 'Confidence (%)';
            yl = [40 100];
            yl2 = [50 100];
            f4 = categorical(squeeze(exp_mat_sess(1:90,1,:,1)));

        case 3
            X       = 100*conf_mat(1:90,:,:)-100*corr_mat(1:90,:,:);
            yTTL = 'Overconfidence (%)';
            yl = [-20 30];
            yl2 = [-50 75];
            f4 = categorical(squeeze(exp_mat_sess(1:90,1,:,1)));
    end
    
    %%% stats %%%
    Xx = squeeze(nanmean(X(:,:,1:24),3));

    n_sub = size(Xx,1);
    vnames = {'GainLoss','Completepartial','subject'};
    f1 = [zeros(2*n_sub,1);ones(2*n_sub,1)];
    f2 = repmat([zeros(n_sub,1);ones(n_sub,1)],2,1);
    rfx = (1:n_sub)';
    f3 = repmat(rfx,4,1);
    
    % glme %
    varNames = {'y','Valence','Info','Participant','Experiment'};
    Xx = Xx./100;
    Xx(Xx==1)=0.9999;
    Xx(Xx==0)= 0.0001;
    idc = find(~isnan(Xx(:)));
    tbl = table(Xx(idc),f1(idc),f2(idc),f3(idc),f4(idc),'VariableNames',varNames);
    
   
    thisReg = fitglme(tbl,'y~ Valence * Info + 1 + Experiment + (1|Participant)','Link','logit');
    thisRegLin = fitglme(tbl,'y~ Valence * Info + 1 + Experiment + (1|Participant)');
    %%% store regression results %%%
    regs{k_meas} = thisReg;
    regsLin{k_meas} = thisRegLin;

    % anova %
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
        
        pExp = floor(log10(pF));
        pDigits = pF./10^pExp;
        if pF<.001
            sigSymb = '***';
        elseif pF<.01
            sigSymb = '**';
        elseif pF<.05
            sigSymb = '*';
        else
            sigSymb = '';
        end
        if pExp<-2
            pStr = [num2str(round(pDigits,2)),'e',num2str(pExp)];
            
            %                    pStr = [num2str(round(pDigits,2)),'x10^',num2str(pExp)];
        else
            pStr = num2str(round(pF,3));
        end
        Ftable{k_meas,iF} =['F = ',num2str(round(F,2)), ', p =',pStr,sigSymb];
        
    end
    
    
    %%% plot partial conditions
    subplot(3,3,sub2ind([3,3],1,k_meas))
    hold on
    for k = [1 3]
        mtp = squeeze(nanmean(X(:,k,:),1));
        stp = squeeze(nanstd(X(:,k,:),0,1))./sqrt(size(X,1));
        fill(Xfill,[mtp'-stp',fliplr(mtp'+stp')],.65*col_mat(k,:),'EdgeColor','none')
        alpha(0.5)
        plot(mtp,'Color',col_mat(k,:));
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
        mtp = squeeze(nanmean(X(:,k,:),1));
        stp = squeeze(nanstd(X(:,k,:),0,1))./sqrt(size(X,1));
        fill(Xfill,[mtp'-stp',fliplr(mtp'+stp')],.65*col_mat(k,:),'EdgeColor','none')
        alpha(0.5)
        plot(mtp,'Color',col_mat(k,:));
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
    dum = squeeze(nanmean(Xx,3))*100;
    dotproperties.useJitter = 1;
    dotproperties.plotLines = [1,2;3,4];
    pirateplot(dum(:,[1,3,2,4])',col_mat([1,3,2,4],:),yl(1),yl(2),12,'','', yTTL,dotproperties);
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

saveas(gcf,'Plots/behaviorLearn.svg');

rmpath('helperfuncs')


%% t-test avg overconf
overconf = nanmean(conf_mat,3)-nanmean(corr_mat,3)*100;
overconfSubAvg = nanmean(overconf,2);
[~,pT,~,t] = ttest(overconfSubAvg);

[~,pTRP,~,tRP] = ttest(overconf(:,1));

%% write anova tables 3
measnames = {'Accuracy','Confidence','Overconfidence'};
dumtable = cell2table(Ftable,'VariableNames', {'Valence','Information','Interaction'},'RowNames',measnames); 
writetable(dumtable,['Results\ANOVAsDataLearning','.csv'],'WriteRowNames',true);
