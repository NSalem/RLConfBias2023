function [h1,h2] = makePlotSsims(outSorted,outConfLSorted,prefMat,accMat,confChoiceMat)

col_mat = [0,1,0;...
    0,1,1;...
    1,0,0;...
    1,0,1];


Xfill = [1:24,fliplr(1:24)];

h1 = figure()
subplot(3,3,1)
hold on
% X = abs(dQ_mat);

X = squeeze(nanmean(outSorted.dQ,3));

for k = [1:4]
    mtp = squeeze(nanmean(X(:,k,:),1));
    stp = squeeze(nanstd(X(:,k,:),0,1))./sqrt(size(X,1));
    fill(Xfill,[mtp'-stp',fliplr(mtp'+stp')],.65*col_mat(k,:),'EdgeColor','none')
    alpha(0.5)
    plot(mtp,'Color',col_mat(k,:));
    xlim([1,24])
    ylim([0,1])
end
ylabel('\DeltaQ')


subplot(3,3,2)
hold on
X = squeeze(nanmean(outSorted.pc,3));
for k = [1:4]
    mtp = squeeze(nanmean(X(:,k,:),1));
    stp = squeeze(nanstd(X(:,k,:),0,1))./sqrt(size(X,1));
    fill(Xfill,[mtp'-stp',fliplr(mtp'+stp')],.65*col_mat(k,:),'EdgeColor','none')
    alpha(0.5)
    plot(mtp,'Color',col_mat(k,:));
    xlim([1,24])
    ylim([0.5,1])
end
ylabel('p(correct)')

subplot(3,3,4)
hold on
plot([1,24],[0,0],'LineStyle','--','Color',[.5,.5,.5,0.5])
X = squeeze(nanmean(outSorted.Qc,3));
for k = [1:4]
    mtp = squeeze(nanmean(X(:,k,:),1));
    stp = squeeze(nanstd(X(:,k,:),0,1))./sqrt(size(X,1));
    fill(Xfill,[mtp'-stp',fliplr(mtp'+stp')],.65*col_mat(k,:),'EdgeColor','none')
    alpha(0.5)
    plot(mtp,'Color',col_mat(k,:));
    xlim([1,24])
    ylim([-1,1])
end
ylabel('Q_C')

subplot(3,3,5)
hold on
plot([1,24],[0,0],'LineStyle','--','Color',[.5,.5,.5,0.5])
X = squeeze(nanmean(outSorted.sigmaQ,3));
for k = [1:4]
    mtp = squeeze(nanmean(X(:,k,:),1));
    stp = squeeze(nanstd(X(:,k,:),0,1))./sqrt(size(X,1));
    fill(Xfill,[mtp'-stp',fliplr(mtp'+stp')],.65*col_mat(k,:),'EdgeColor','none')
    alpha(0.5)
    plot(mtp,'Color',col_mat(k,:));
    xlim([1,24])
    ylim([-1,1])
    
end
ylabel('\SigmaQ')

subplot(3,3,6)
hold on
plot([1,24],[0,0],'LineStyle','--','Color',[.5,.5,.5,0.5])
X = squeeze(nanmean(outSorted.V,3));
for k = [1:4]
    mtp = squeeze(nanmean(X(:,k,:),1));
    stp = squeeze(nanstd(X(:,k,:),0,1))./sqrt(size(X,1));
    fill(Xfill,[mtp'-stp',fliplr(mtp'+stp')],.65*col_mat(k,:),'EdgeColor','none')
    alpha(0.5)
    plot(mtp,'Color',col_mat(k,:));
    xlim([1,24])
    ylim([-1,1])
    
end
ylabel('V')

idmodels = [2:4];
for iconfmodelLT = 1:numel(idmodels)
    switch iconfmodelLT
        case 1
            yL = 'Confidence Q_C model'
        case 2
            yL = 'Confidence \SigmaQ model'
        case 3
            yL = 'Confidence V model'
    end
    subplot(3,3,iconfmodelLT+6)
    hold on
    X = squeeze(nanmean(nanmean(outConfLSorted(:,:,idmodels(iconfmodelLT),:,:,:),4),1));
    for k = [1:4]
        mtp = squeeze(nanmean(X(:,k,:),1));
        stp = squeeze(nanstd(X(:,k,:),0,1))./sqrt(size(X,1));
        fill(Xfill,[mtp'-stp',fliplr(mtp'+stp')],.65*col_mat(k,:),'EdgeColor','none')
        alpha(0.5)
        plot(mtp,'Color',col_mat(k,:));
%         ylim([60,90])
        ylim([.5,1])
        xlim([1,24])
    end
    ylabel(yL)
    xlabel('Trial')
end
    
    
    
%% Transfer
% confChoiceMat = 
%% post test matrices
h2 = figure()
for iplot = 1:6
    switch iplot
        case 1
            myMat = squeeze(nanmean(prefMat(:,:,:),1))'/4*100;
            clim = [0,100];
            T = 'p(choose S1)';
        case 2
            myMat = squeeze(nanmean(accMat(:,:,:),1))'/4*100;
            clim = [0,100];
            T = 'Accuracy (%)';
        case 3 
            continue
        case 4
            myMat = squeeze(nanmean(confChoiceMat(2,:,:,:),2))'*100;
            clim = [50,95];
            T = 'Confidence  (Q_c model)';
        case 5
            myMat = squeeze(nanmean(confChoiceMat(3,:,:,:),2))'*100;
            clim = [50,95];
            T = 'Confidence (\SigmaQ model)';
        case 6
            myMat = squeeze(nanmean(confChoiceMat(4,:,:,:),2))'*100;
            clim = [50,95];
            T = 'Confidence (V model)';
    end
        subplot(2,3,iplot)

    if iplot <3 %% for preference and accuracy symmetric matrix. For (over)confidence we distinguish chosen/unchosen
        myMat = tril(myMat);
    end
    
    myMat(myMat==0) = NaN;
    b = imagesc(myMat);
    set(b,'AlphaData',~isnan(myMat));
    
    l = colorbar;
%     caxis(clim)
    ylabel(l,T,'FontSize',11)
    
    newIdc = [6,5,2,1,8,7,4,3];
    lg={'\color{green}G_{75}','\color{green}G_{25}','\color{cyan}G_{75}','\color{cyan}G_{25}',...
        '\color{red}L_{25}','\color{red}L_{75}','\color{magenta}L_{25}','\color{magenta}L_{75}'};
    lg={lg{newIdc}};
    xticks(1:8)
    yticks(1:8)
    xticklabels(lg)
    yticklabels(lg)
    set(gca,'FontSize',8)
    
    if iplot>2
        xlabel('Unchosen')
        ylabel('Chosen')
    else
        xlabel('S1')
        ylabel('S2')
    end
end

set(gcf,'Position',[0,0,1200,600])

