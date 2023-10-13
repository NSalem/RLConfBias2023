%Plots to evaluate  model identifiability and parameter
%recovery for choice, using results from modelIdentRecovChoice.m
%(Palminteri et al 2016,
%https://doi.org/10.1016/j.tics.2017.03.011)
%Plots adapted from Correa et al 2018
%(https://doi.org/10.1523/JNEUROSCI.0457-18.2018)

% load('Results\SIMS_ResultsIdentRecov.mat')


% files ={'Results\SIMS_ResultsIdentRecov50c.mat','Results\SIMS_ResultsIdentRecov.mat'};
% outnames = {'50c','10c100c'};

files ={'Results\ResultsIdentRecovChoice.mat',};
outnames = {'10c100c_2022'};

for ifile = 1
    
    load(files{ifile})
    
%     modellabels = {'ABS','REL_0','REL_{Q_U}','REL_{OTHER}','REL_{LAST}',...
%         'CONF','CONFREL_0','CONFREL_{Q_U}','CONFREL_{OTHER}','CONFREL_{LAST}'};
%     
        modellabels = {'ABS','REL_0','REL_{Q_U}','REL_{OTHER}','REL_{LAST}','REL_{V}','REL_{Q+V}'...
        'CONF','CONFREL_0','CONFREL_{Q_U}','CONFREL_{OTHER}','CONFREL_{LAST}','CONFREL_{V}','CONFREL_{Q+V}'};
    
    
    nparams = size(Rest,1);

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
                whichmodels = [1:5,8:12];
        end
    nmodels = size(pxpMat,1);
    nsims = size(pxpMat,3);
    mean_pxp     = squeeze(mean(pxpMat,3));
    n_goodcl    = bmMat;
    
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
        set(gca,'XTick',1:nmodels,...
            'YTick',1:nmodels,...
            'XTickLabel',{modellabels{whichmodels}},...
            'YTickLabel',{modellabels{fliplr(whichmodels)}},...
            'FontName','Arial');
        ylabel('Simulated model')
        xlabel('Estimated model')
        xtickangle(90)
        c = colorbar;
        c.Label.String = lbl;
    end
    
    saveas(gcf,['Plots/identRecovFig1_',outnames{ifile},'_',label,'.svg'])
    end
    %% Fig 2
    h2 = figure('Units', 'pixels', ...
        'Position', [400 150 600 600]);
    set(h2,'Color',[1,1,1])
    
    LAB = {'\beta','\alpha_{CON}','\alpha_{DIS}','\alpha_V','w'};
    
    for k = 1:5
        
        subplot(3,2,k)
        hold on
        if  k ==1
            x = 0:0.2:10;
            distr_tp = gampdf(x,1.25,5);
            xl = [0 10];
        else
            x = 0:0.01:1;
            distr_tp = betapdf(x,1.1,1.1);
            xl = [0 1];
        end
        
        ax1 = gca; % current axes
        
        plot(ax1,x,distr_tp,'Color',.5*[1,1,1])
        set(ax1,'XLim',xl,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'XTickLabel',[],...
            'YTickLabel',[],...
            'FontName','Arial')
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
        
        
        X = squeeze(genpar(:,k,:));
        Y = squeeze(recpar(:,k,:));
        
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
    
    saveas(gcf,['Plots/identRecovFig2_',outnames{ifile},'.svg'])

    
    %% Fig 3
    figure('Units', 'pixels', ...
        'Position', [400 200 200 350],...
        'Color',[1,1,1]);
    
    mtp = NaN(2,nparams);
    stp = NaN(2,nparams);
    for k = 1:nparams
        
        
        mtp = STORE_REG(k).stats.beta;
        stp = STORE_REG(k).stats.se;
        subplot(3,2,k)
        hold on
        bar(mtp,'FaceColor',.9.*[1,1,1])
        errorbar(mtp,stp,'k','LineStyle','none')
        set(gca,'YLim',[-0.5 1.5],...
            'XLim',[0 3],...
            'XTick',[1 2],...
            'XTickLabel',{'\beta_0','\beta_1'},'FontName','Arial')
        ylabel(LAB{k});
        
    end

    saveas(gcf,['Plots/identRecovFig3_',outnames{ifile},'.svg'])

    %% Fig 4
    for k = 1:2
        
        figure('Units', 'pixels', ...
            'Position', [400 200 450 350]);
        set(gcf,'Color',[1,1,1])
        switch k
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
        ylabel('Parameter')
        xlabel('Parameter')
        set(gca,'XTickLabel',LAB,'XTick',[1:size(mat_tp,1)],...
            'YTick',[1:size(mat_tp,1)],...
            'YTickLabel',fliplr(LAB),'FontName','Arial')
        c = colorbar;
        c.Label.String = lbl;
        caxis(cax)
        
        if k == 1
            for ir = 1:size(mat_tp,1)
                for ic = 1:size(mat_tp,2)
                    if (size(mat_tp)-ir+1)<=ic
                        text(ic,ir, sprintf('%0.2f',(mat_tp(ir,ic))), 'FontSize', 11,'HorizontalAlignment','center')
                    end
                end
            end
        end
        
    end
    
    rDiag = eye(nparams).*rmat_avg;
    rDiag(rDiag==0) = NaN;
    nanmean(rDiag(:))
    
    rNdiag = (ones(nparams)-eye(nparams)).*rmat_avg;
    rNdiag(rNdiag==0) = NaN;
    nanmean(rNdiag(:))
    
    r2Diag = eye(nparams).*r2mat_avg;
    r2Diag(r2Diag==0) = NaN;
    nanmean(r2Diag(:))
    
    r2Ndiag = (ones(nparams)-eye(nparams)).*r2mat_avg;
    r2Ndiag(r2Ndiag==0) = NaN;
    nanmean(r2Ndiag(:))
    5
    if k ==1
    saveas(gcf,['Plots/identRecovFig4_',outnames{ifile},'.svg'])
    elseif k==2
        saveas(gcf,['Plots/identRecovFig5_',outnames{ifile},'.svg'])    
    end
end