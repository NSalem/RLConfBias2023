%Plots to evaluate  model identifiability and parameter
%recovery for confidence, using results from modelIdentRecovConfidence.m
%(Palminteri et al 2016,
%https://doi.org/10.1016/j.tics.2017.03.011)
%Plots adapted from Correa et al 2018
%(https://doi.org/10.1523/JNEUROSCI.0457-18.2018)

load('Results\ResultsIdentRecovConfidence.mat')
modelnamesLT = {'No bias','Q_C','\SigmaQ','V','Q_C+V','\SigmaQ+V','conf_{t-1}','Q_C+conf_{t-1}','\SigmaQ+conf_{t-1}','V+conf_{t-1}','Q_C+V+conf_{t-1}','\SigmaQ+V+conf_{t-1}','Q_C+Q_U+V','Q_C+Q_U+V+conf_{t-1}'};
modelnamesTT = {'No bias','Q_C','\SigmaQ','V','Q_C+V','\SigmaQ+V','conf_{t-1}','Q_C+conf_{t-1}','\SigmaQ+conf_{t-1}','V+conf_{t-1}','Q_C+V+conf_{t-1}','\SigmaQ+V+conf_{t-1}','Q_C+Q_U+V','Q_C+Q_U+V+conf_{t-1}'};

nsims = size(pxpLT,3);
nmodelsLT = size(pxpLT,1);
nparamsLT = size(RestLT,1);
nmodelsTT = size(pxpTT,1);
nparamsTT = size(RestTT,1);

%% plot
for itask = 1:2 %learning and transfer
    switch itask
        case 1
            bm = bmLT;
            bmReduced = bmLTReduced;
            pxp = pxpLT;
            pxpReduced = pxpLTReduced;
            genpar = genparLT;
            recpar = recparLT;
%             nmodels = nmodelsLT;
            modelnames = modelnamesLT;
            Rest = RestLT;
            R2est = R2estLT;
            Rest2 = RestLT2;
            R2est2 = R2estLT2;
            npars = nparamsLT;
            genpar2 = genparLT2;
            recpar2 = recparLT2;
            npars2 = size(genparLT2,2);
            LAB = {'B_0','B_{|\DeltaQ|}','B_{Q_c}','B_{Conf(t-1)}'};
            LAB2 = {'B_0','B_{Q_C}','B_{Q_U}','B_{V}','B_{Conf(t-1)}'};
            
        case 2
            bm = bmTT;
            pxp = pxpTT;
            bmReduced = bmTTReduced;
            pxpReduced = pxpTTReduced;
            genpar = genparTT;
            recpar = recparTT;
%             nmodels = nmodelsTT;
            modelnames = modelnamesTT;
            Rest = RestTT;
            R2est = R2estTT;
            Rest2 = RestTT2;
            R2est2 = R2estTT2;
            npars = nparamsTT;
            
            genpar2 = genparTT2;
            recpar2 = recparTT2;
            npars2 = size(genparTT2,2);
            
            LAB = {'B_0','B_{|\DeltaQ|}','B_{Q_c}'};
            LAB2 = {'B_0','B_{Q_C}','B_{Q_U}','B_{V}'};
            
    end
    
    

    for modelsel = 1:2
        switch modelsel
            case 1 
                pxpMat = pxp;
                bmMat = bm;
                label = 'full';
                whichmodels = [1:14];
            case 2
                pxpMat = pxpReduced;
                bmMat = bmReduced;
                label = 'reduced'
                if itask ==1
                whichmodels = [7:10];
                else 
                    whichmodels = 1:4;
                end
        end
        
    % compute confusion matrices
    mean_pxp     = squeeze(mean(pxpMat,3));
    n_goodcl    = squeeze(sum(bmMat,3));
    
    nmodels = size(pxpMat,2);
    
    % Fig 1
    h1 = figure('Units', 'pixels', ...
        'Position', [400 200 1000 350]);
    set(h1,'Color',[1,1,1])
    
    for k = 1:2
        subplot(1,2,k)
        
        switch k
            case 1
                mtp = mean_pxp;
                lbl = 'Exceedance probability (%)';
            case 2
                mtp = n_goodcl/nsims*100;
                lbl = '% Best model';
        end
        
        colormap(flipud(gray))
        imagesc(flipud(mtp))
        ylabel('Simulated model')
        xlabel('Estimated model')
        set(gca,'XTick',1:nmodels,...
            'YTick',1:nmodels,...
            'XTickLabel',{modelnames{whichmodels}},...
            'YTickLabel',fliplr({modelnames{whichmodels}}))
        xtickangle(90)
        c = colorbar;
        c.Label.String = lbl;
    end
    
    saveas(gcf,['Plots/identRecovConfFig1_',label,'.svg'])
    end
    %% Fig 2
    h2 = figure('Units', 'pixels', ...
        'Position', [400 150 600 600]);
    set(h2,'Color',[1,1,1])
    
    for k = 1:npars
        X = squeeze(genpar(:,k,:));
        Y = squeeze(recpar(:,k,:));
        xl = [min([X(:);Y(:)]),max([X(:);Y(:)])];
        
        subplot(3,2,k)
        hold on
        ax1 = gca; % current axes
        set(ax1,'XLim',xl,...
            'Ylim',xl,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'XTickLabel',[],...
            'YTickLabel',[],'FontName','Arial')
        ax1_pos = ax1.Position;
        ax2 = axes('Position',ax1_pos,...
            'YLim',xl,...
            'XAxisLocation','bottom',...
            'YAxisLocation','left',...
            'Color','none',...
            'XLim',xl);
        set(ax2,'FontName','Arial')
        xlabel(strcat(['Simulated ',LAB{k}]));
        ylabel(strcat(['Estimated ',LAB{k}]));
        
        hold on
        
        plot(ax2,xl,xl,':k',...
            'LineWidth',2)
        plot(ax2,X(:),Y(:),'o',...
            'MarkerFaceColor',[1,1,1],...
            'MarkerEdgeColor',[0,0,0])
        
        [b,~,stats]    = glmfit(X(:),Y(:),'normal');
        STORE_REG(k).b = b;
        STORE_REG(k).stats = stats;
        
        XX             = linspace(min(X(:)),max(X(:)),1000);
        [Yf,Yl,Yh]      = glmval(b,XX,'identity',stats,'confidence',0.95);
        XXX            = sortrows([XX',Yf,Yf-Yl,Yf+Yh],1);
        
        Xfill = [XXX(:,1);flipud(XXX(:,1))];
        fill(Xfill,[XXX(:,3);flipud(XXX(:,4))],.7*[1,1,1],'EdgeColor','none',...
            'Parent', ax2)
        alpha(0.5)
        hModel         =   plot(ax2,XXX(:,1),XXX(:,2),'-',...
            'Color',.5*[1,0,0],...
            'LineWidth',2)
        set(ax2,'Position',ax1_pos)
    end
    
    saveas(gcf,['Plots/identRecovConfFig2.svg'])
    
    
    %% Fig 3
    figure('Units', 'pixels', ...
        'Position', [400 200 600 200],...
        'Color',[1,1,1]);
    
    mtp = NaN(2,npars);
    stp = NaN(2,npars);
    for k = 1:npars
        mtp = STORE_REG(k).stats.beta;
        stp = STORE_REG(k).stats.se;
        subplot(1,5,k)
        hold on
        bar(mtp,'FaceColor',.9.*[1,1,1])
        errorbar(mtp,stp,'k','LineStyle','none')
        set(gca,'YLim',[-0.5 2],...
            'XLim',[0 3],...
            'XTick',[1 2],...
            'XTickLabel',{'intercept','slope'},...
            'FontName','Arial')
        title(LAB{k});
        xtickangle(90);
    end
    saveas(gcf,['Plots/identRecovConfFig3.svg'])
    
    %% Fig 4
    for j = 1:2
        figure('Units', 'pixels', ...
            'Position', [400 200 450 350]);
        set(gcf,'Color',[1,1,1])
        switch j
            case 1
                colormap(parula)
                mat_tp = squeeze(mean(Rest,3));
                lbl = 'Pearson Correlation (R)';
                cax = [-1 1];
                rmat_avg = mat_tp;
            case 2
                colormap(flipud(gray))
                mat_tp = squeeze(mean(R2est,3));
                lbl = 'Explained variance (R^2)';
                cax = [0 1];
                r2mat_avg = mat_tp;
        end
        
        mat_tp = flipud(mat_tp);
        imagesc(mat_tp)
        set(gca,'XTickLabel',LAB,...
            'XTick',1:numel(LAB),...
            'YTick',1:numel(LAB),...
            'YTickLabel',fliplr(LAB),'FontName','Arial')
        ylabel('Parameter')
        xlabel('Parameter')
        c = colorbar;
        c.Label.String = lbl;
        caxis(cax)
        
        if j == 1
            for ir = 1:size(mat_tp,1)
                for ic = 1:size(mat_tp,2)
                    if (size(mat_tp)-ir+1)<=ic
                        text(ic,ir, sprintf('%0.2f',(mat_tp(ir,ic))), 'FontSize', 11,'HorizontalAlignment','center','FontName','Arial')
                    end
                end
            end
        end
        
        if j ==1
            saveas(gcf,['Plots/identRecovConfFig4.svg'])
        elseif j==2
            saveas(gcf,['Plots/identRecovConfFig5.svg'])
        end
        
    end
    
    
    %%% figures for Qc+Qu+V
    %% Fig 2  Qc+Qu+V
    h2 = figure('Units', 'pixels', ...
        'Position', [400 150 600 600]);
    set(h2,'Color',[1,1,1])
    
    for k = 1:npars2
        X = squeeze(genpar2(:,k,:));
        Y = squeeze(recpar2(:,k,:));
        xl = [min([X(:);Y(:)]),max([X(:);Y(:)])];
        
        subplot(3,2,k)
        hold on
        ax1 = gca; % current axes
        set(ax1,'XLim',xl,...
            'Ylim',xl,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'XTickLabel',[],...
            'YTickLabel',[],'FontName','Arial')
        ax1_pos = ax1.Position;
        ax2 = axes('Position',ax1_pos,...
            'YLim',xl,...
            'XAxisLocation','bottom',...
            'YAxisLocation','left',...
            'Color','none',...
            'XLim',xl);
        set(ax2,'FontName','Arial')
        xlabel(strcat(['Simulated ',LAB2{k}]));
        ylabel(strcat(['Estimated ',LAB2{k}]));
        
        hold on
        
        plot(ax2,xl,xl,':k',...
            'LineWidth',2)
        plot(ax2,X(:),Y(:),'o',...
            'MarkerFaceColor',[1,1,1],...
            'MarkerEdgeColor',[0,0,0])
        
        [b,~,stats]    = glmfit(X(:),Y(:),'normal');
        STORE_REG(k).b = b;
        STORE_REG(k).stats = stats;
        
        XX             = linspace(min(X(:)),max(X(:)),1000);
        [Yf,Yl,Yh]      = glmval(b,XX,'identity',stats,'confidence',0.95);
        XXX            = sortrows([XX',Yf,Yf-Yl,Yf+Yh],1);
        
        Xfill = [XXX(:,1);flipud(XXX(:,1))];
        fill(Xfill,[XXX(:,3);flipud(XXX(:,4))],.7*[1,1,1],'EdgeColor','none',...
            'Parent', ax2)
        alpha(0.5)
        hModel         =   plot(ax2,XXX(:,1),XXX(:,2),'-',...
            'Color',.5*[1,0,0],...
            'LineWidth',2)
        set(ax2,'Position',ax1_pos)
    end
    
    saveas(gcf,['Plots/identRecovConfQcQuVFig2.svg'])
    
    
    %% Fig 3
    figure('Units', 'pixels', ...
        'Position', [400 200 600 200],...
        'Color',[1,1,1]);
    
    mtp = NaN(2,npars2);
    stp = NaN(2,npars2);
    for k = 1:npars2
        mtp = STORE_REG(k).stats.beta;
        stp = STORE_REG(k).stats.se;
        subplot(1,5,k)
        hold on
        bar(mtp,'FaceColor',.9.*[1,1,1])
        errorbar(mtp,stp,'k','LineStyle','none')
        set(gca,'YLim',[-0.5 2],...
            'XLim',[0 3],...
            'XTick',[1 2],...
            'XTickLabel',{'intercept','slope'},...
            'FontName','Arial')
        title(LAB2{k});
        xtickangle(90);
    end
    saveas(gcf,['Plots/identRecovConfQcQuVFig3.svg'])
    
    %% Fig 4 Qc+Qu+V model
    for j = 1:2
        figure('Units', 'pixels', ...
            'Position', [400 200 450 350]);
        set(gcf,'Color',[1,1,1])
        switch j
            case 1
                colormap(parula)
                mat_tp = squeeze(mean(Rest2,3));
                lbl = 'Pearson Correlation (R)';
                cax = [-1 1];
                rmat_avg = mat_tp;
            case 2
                colormap(flipud(gray))
                mat_tp = squeeze(mean(R2est2,3));
                lbl = 'Explained variance (R^2)';
                cax = [0 1];
                r2mat_avg = mat_tp;
        end
        
        mat_tp = flipud(mat_tp);
        imagesc(mat_tp)
        set(gca,'XTickLabel',LAB2,...
            'XTick',1:numel(LAB2),...
            'YTick',1:numel(LAB2),...
            'YTickLabel',fliplr(LAB2),'FontName','Arial')
        ylabel('Parameter')
        xlabel('Parameter')
        c = colorbar;
        c.Label.String = lbl;
        caxis(cax)
        
        if j == 1
            for ir = 1:size(mat_tp,1)
                for ic = 1:size(mat_tp,2)
                    if (size(mat_tp)-ir+1)<=ic
                        text(ic,ir, sprintf('%0.2f',(mat_tp(ir,ic))), 'FontSize', 11,'HorizontalAlignment','center','FontName','Arial')
                    end
                end
            end
        end
        
        if j ==1
            saveas(gcf,['Plots/identRecovConfQcQuVFig4.svg'])
        elseif j==2
            saveas(gcf,['Plots/identRecovConfQcQuVFig5.svg'])
        end
        
    end
    
end


