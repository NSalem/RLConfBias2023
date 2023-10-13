clear

f = load('Results\reg_conflogit_learning.mat','whichLearnModel' ,'isConfPrev' ,'whichConfModel');
whichLearnModel = f.whichLearnModel; isConfPrev = f.isConfPrev; whichConfModel = f.whichConfModel;

idcmodels = 1:numel(whichLearnModel)%find(~modelvars.isConfPrev);
whichmodelsLT = idcmodels;
whichmodelsTT = idcmodels;
nmodelsLT = numel(whichmodelsLT);
nmodelsTT = numel(whichmodelsTT);

% load results if already existent
if exist('Results\SIMS_CONFIDENCE_RsultsIdentRecovALL.mat','file')
    doIdentRecov = false;
    load ('Results\SIMS_CONFIDENCE_RsultsIdentRecovALL.mat')
    nsims = size(bmLT,3);
else
    load('Results\SIMS_conf2022_04_08_CONFIDENCELOGIT.mat');
    nsims = isim;
    doIdentRecov = true;
    bmLT      = zeros(nmodelsLT,nmodelsLT,nsims);  % best model
    pxpLT     = zeros(nmodelsLT,nmodelsLT,nsims);  % Protected exceedance probability
    bmTT      = zeros(nmodelsTT,nmodelsTT,nsims);  % best model
    pxpTT     = zeros(nmodelsTT,nmodelsTT,nsims);  % Protected exceedance probability
end

if doIdentRecov
       
    nsubs = size(regress.coeffs,2);
    
    %% prepare parameters for recovery    
    % find index and number of parameters for best model (with
    % autocorrelation parameter)
    model = find(whichLearnModel(idcmodels) == 11 & whichConfModel(idcmodels) == 2 & isConfPrev(idcmodels));
    
    nparamsLT = numel(genparamsconf{1,1,model}); %number of parameters of full model
    nparamsTT = numel(genparamsconf{1,1,model}); %number of parameters of full model
    
    %% get matrices generative and recovered parameters (sub x par x sim)
    genparLT = nan(nsubs,nparamsLT,nsims);
    recparLT = nan(nsubs,nparamsLT,nsims);
    genparTT = nan(nsubs,nparamsLT,nsims);
    recparTT = nan(nsubs,nparamsLT,nsims);
    
    for isub  = 1:nsubs
        for isim = 1:nsims
            genparLT(isub,:,isim) = genparamsconf{isim,isub,model};
            recparLT(isub,:,isim) = [regress.coeffs{isim,isub,model,model};regress.RMSE(isim,isub,model,model)];
            genparTT(isub,:,isim) = genparamsconfpost{isim,isub,model};
            recparTT(isub,:,isim) = [regress.coeffsPost{isim,isub,model,model};regress.RMSEPost(isim,isub,model,model)];
        end
    end
    
    RestLT    = NaN(nparamsLT,nparamsLT,nsims);
    R2estLT   = NaN(nparamsLT,nparamsLT,nsims);
    RestTT    = NaN(nparamsTT,nparamsTT,nsims);
    R2estTT   = NaN(nparamsTT,nparamsTT,nsims);
    
    %% do BMC and parameter recovery matrices
    options = struct();
    criterionLT = regress.BIC;
    criterionTT = regress.BICPost;
    options.DisplayWin = 0;
    
    for isim = 1:nsims
        for igenmodel = 1:numel(whichmodelsLT)
            genmodel = whichmodelsLT(igenmodel);
            [~, out] = VBA_groupBMC(squeeze(-criterionLT(isim,:,genmodel,whichmodelsLT))'/2,options);
            [~,winningmodel] = max(out.pxp);
            pxpLT(igenmodel,:,isim) = out.pxp;
            bmLT(igenmodel,winningmodel,isim) = 1;
        end
        for igenmodel = 1:numel(whichmodelsTT)
            genmodel = whichmodelsTT(igenmodel);
            
            [~, out] = VBA_groupBMC(squeeze(-criterionTT(isim,:,genmodel,whichmodelsTT))'/2,options);
            [~,winningmodel] = max(out.pxp);
            pxpTT(igenmodel,:,isim) = out.pxp;
            bmTT(igenmodel,winningmodel,isim) = 1;
        end   
        
        % compute correlations between parameters used to simulate the data,
        % and recovered (i.e. estimated) parameters
        RestLT(:,:,isim) = corr(squeeze(recparLT(:,:,isim)));
        RestTT(:,:,isim) = corr(squeeze(recparTT(:,:,isim)));
        for ipar = 1:nparamsLT
            RestLT(ipar,ipar,isim) = corr(squeeze(recparLT(:,ipar,isim)),squeeze(genparLT(:,ipar,isim)));
            RestTT(ipar,ipar,isim) = corr(squeeze(recparTT(:,ipar,isim)),squeeze(genparTT(:,ipar,isim)));
        end
        R2estLT(:,:,isim) = RestLT(:,:,isim).*RestLT(:,:,isim);
        R2estTT(:,:,isim) = RestTT(:,:,isim).*RestTT(:,:,isim);
    end
    
    %%% save
    infostr = struct('nsim','nmodels','idcmodels','model')
    save('Results\SIMS_CONFIDENCE_RsultsIdentRecovALL.mat','bmLT','bmTT',...
        'pxpLT','pxpTT','RestLT','RestTT','R2estLT','R2estTT',...
        'genparLT','genparTT','recparLT','recparTT','infostr')
    
end

%% plot
for k = 1:2 %learning and transfer
    switch k
        case 1
            bm = bmLT;
            pxp = pxpLT;
            genpar = genparLT;
            recpar = recparLT;
            nmodels = nmodelsLT;
            Rest = RestLT;
            R2est = R2estLT;
        case 2
            bm = bmTT;
            pxp = pxpTT;
            genpar = genparLT;
            recpar = recparLT;
            nmodels = nmodelsTT;
            R2est = R2estTT;
            R2est = R2estTT;
    end
    
    
    % compute confusion matrices
    mean_pxp     = squeeze(mean(pxp,3));
    n_goodcl    = squeeze(sum(bm,3));
    
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
                mtp = n_goodcl;
                lbl = '% Best model';
        end
        
        colormap(flipud(gray))
        imagesc(flipud(mtp))
        ylabel('simulated model #')
        xlabel('estimated model #')
        set(gca,'XTick',1:nmodels,...
            'YTick',1:nmodels,...
            'XTickLabel',(1:nmodels),...
            'YTickLabel',fliplr(1:nmodels))
        c = colorbar;
        c.Label.String = lbl;
    end
   
    %% Fig 2
    h2 = figure('Units', 'pixels', ...
        'Position', [400 150 600 600]);
    set(h2,'Color',[1,1,1])
    
    LAB = {'b0','b_|\DeltaQ|','b_{Q_c}','b_{Conf(t-1)}','RMSE'};
    
    for k = 1:nparamsLT
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
            'YTickLabel',[])
        ax1_pos = ax1.Position;
        ax2 = axes('Position',ax1_pos,...
            'YLim',xl,...
            'XAxisLocation','bottom',...
            'YAxisLocation','left',...
            'Color','none',...
            'XLim',xl);
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
    
    
    
    %% Fig 3
    figure('Units', 'pixels', ...
        'Position', [400 200 200 350],...
        'Color',[1,1,1]);
    
    mtp = NaN(2,nparamsLT);
    stp = NaN(2,nparamsLT);
    for k = 1:nparamsLT
        
        
        mtp = STORE_REG(k).stats.beta;
        stp = STORE_REG(k).stats.se;
        subplot(3,2,k)
        hold on
        bar(mtp,'FaceColor',.9.*[1,1,1])
        errorbar(mtp,stp,'k','LineStyle','none')
        set(gca,'YLim',[-0.5 1.5],...
            'XLim',[0 3],...
            'XTick',[1 2],...
            'XTickLabel',{'\beta_0','\beta_1'})
        ylabel(LAB{k});
        
    end
    
    
    
    
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
        ylabel('Parameter')
        xlabel('Parameter')
        set(gca,'XTickLabel',LAB,...
            'YTickLabel',fliplr(LAB))
        c = colorbar;
        c.Label.String = lbl;
        caxis(cax)
        
        if j == 1
            for ir = 1:size(mat_tp,1)
                for ic = 1:size(mat_tp,2)
                    if (size(mat_tp)-ir+1)<=ic
                        text(ic,ir, sprintf('%0.2f',(mat_tp(ir,ic))), 'FontSize', 11,'HorizontalAlignment','center')
                    end
                end
            end
        end
    end
end


