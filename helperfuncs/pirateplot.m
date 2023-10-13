function [Nbar,Nsub,jitter] = pirateplot(DataCell,Colors,Yinf,Ysup,Font,Title,LabelX,LabelY,dotproperties,varargin)

% Sophie Bavard - December 2018
% Creates a violin plot with mean, error bars, confidence interval, kernel density.
% Warning: the function can accept any number of arguments > 9.
% After the Title, LabelX, LabelY : varargin for bar names under X-axis

if ~exist('dotproperties')||~isfield(dotproperties,'alpha')
    dotproperties.alpha = .5;
end
if ~exist('dotproperties') || ~isfield(dotproperties,'useJitter')
    dotproperties.useJitter = 1;
end

if ~exist('dotproperties') || ~isfield(dotproperties,'plotLines')
    dotproperties.plotLines = [];
end
% transforms the Data matrix into cell format if needed
if iscell(DataCell)==0
    DataCell = num2cell(DataCell', 1)';
end

% number of factors/groups/conditions
Nbar = size(DataCell,1);
% bar size
Wbar = 0.75;

% confidence interval
ConfInter = 0.95;

% color of the box + error bar
trace = [0.5 0.5 0.5];


% check number of elements per cell, if they're the same use same jitter

cellSizes = [];
for icell = 1:size(DataCell)
    cellSizes(icell) = numel(DataCell{icell});
end
hold on
if numel(unique(cellSizes))==1
    myMat = cell2mat(DataCell');
    if dotproperties.useJitter
        jitter = .75*(rand(1,cellSizes(1))-0.5);
    else
        jitter = 0*(rand(1,cellSizes(1))-0.5);
    end
end


for n = 1:Nbar
    % calculate kernel density estimation for the violin
    clear DataMatrix
    DataMatrix = DataCell{n,:}';
    Nsub = length(DataMatrix(~isnan(DataMatrix)));
    
    [density, value] = ksdensity(DataMatrix, 'Bandwidth', 0.9 * min(std(DataMatrix), iqr(DataMatrix)/1.34) * Nsub^(-1/5)); % change Bandwidth for violin shape. Default MATLAB: std(DataMatrix)*(4/(3*Nsub))^(1/5)
    density = density(value >= min(DataMatrix) & value <= max(DataMatrix));
    value = value(value >= min(DataMatrix) & value <= max(DataMatrix));
    value(1) = min(DataMatrix);
    value(end) = max(DataMatrix);
    
    % all data is identical
    if min(DataMatrix) == max(DataMatrix)
        density = 1;
    end
    width = Wbar/2/max(density);
    
    % INDIVIDUAL DOTS INSIDE VIOLIN
    if length(density) > 1
        jitterstrength = interp1(value, density*width, DataMatrix);
    else % all data is identical
        jitterstrength = density*width;
    end
    if numel(unique(cellSizes))==1
        jitterstrengthAll(n,:) = jitterstrength;
    end
end


%% plot lines connecting dots
if ~(isempty(dotproperties.plotLines))
    pairMat = dotproperties.plotLines;
    for ipair = 1:size(pairMat,1)
        for isub = 1:size(myMat,1)
            plot(pairMat(ipair,:)+jitter(isub).*jitterstrengthAll(pairMat(ipair,:),isub)',myMat(isub,pairMat(ipair,:)),'Color',[.8,.8,.8]);
        end
    end
end


for n = 1:Nbar
    
    clear DataMatrix
    DataMatrix = DataCell{n,:}';
    
    % number of subjects
    Nsub = length(DataMatrix(~isnan(DataMatrix)));
    
    curve = nanmean(DataMatrix,2);
    sem   = nanstd(DataMatrix')'/sqrt(Nsub);
    conf  = tinv(1 - 0.5*(1-ConfInter),Nsub);
    
    if dotproperties.useJitter && numel(unique(cellSizes))>1
        jitter = 2*(rand(size(DataMatrix))-0.5);
    end
    % error bar rectangle
    % rectangle('Position',[n- Wbar/4, curve-sem, Wbar/2, sem*2],'EdgeColor',[0.5 0.5 0.5]);
    % hold on
    
    % PLOT THE VIOLINS
    
    % calculate kernel density estimation for the violin
    [density, value] = ksdensity(DataMatrix, 'Bandwidth', 0.9 * min(std(DataMatrix), iqr(DataMatrix)/1.34) * Nsub^(-1/5)); % change Bandwidth for violin shape. Default MATLAB: std(DataMatrix)*(4/(3*Nsub))^(1/5)
    density = density(value >= min(DataMatrix) & value <= max(DataMatrix));
    value = value(value >= min(DataMatrix) & value <= max(DataMatrix));
    value(1) = min(DataMatrix);
    value(end) = max(DataMatrix);
    
    % all data is identical
    if min(DataMatrix) == max(DataMatrix)
        density = 1;
    end
    width = Wbar/2/max(density);
    
    % plot the violin
    fill([n+density*width n-density(end:-1:1)*width],...
        [value value(end:-1:1)],...
        Colors(n,:).*0.65,...
        'EdgeColor', 'none',...
        'FaceAlpha',0.3);
    hold on
    
    % INDIVIDUAL DOTS INSIDE VIOLIN
    if length(density) > 1
        jitterstrength = interp1(value, density*width, DataMatrix);
    else % all data is identical
        jitterstrength = density*width;
    end
    %     jitterstrengthAll{:} = jitterstrength;
    %     jitter = 2*(rand(size(DataMatrix))-0.5);
    scatter(n + jitter.*jitterstrength.*dotproperties.useJitter, DataMatrix, 10,...
        'marker','o',...
        'LineWidth',.5,...
        'MarkerFaceColor',Colors(n,:),...
        'MarkerEdgeColor',[1,1,1],...
        'MarkerEdgeAlpha',dotproperties.alpha);
    
    % COLORED BARS
    % bar(n,curve,...
    %     'FaceColor',Colors(n,:),...
    %     'EdgeColor','none',...
    %     'BarWidth',Wbar,...
    %     'LineWidth',1,...
    %     'FaceAlpha',0.15);
    % hold on
    
    % MEAN HORIZONTAL BAR
    xMean = [n - Wbar/2; n + Wbar/2];
    yMean = [curve; curve];
    hold on, plot(xMean,yMean,'-','LineWidth',2,'Color',Colors(n,:).*0.65);
    
    % ERROR BARS
    errorbar(n,curve,sem,...
        'Color',Colors(n,:).*0.8,...%trace-0.1,...
        'LineStyle','none',...
        'LineWidth',1);
    hold on
    
    % CONFIDENCE INTERVAL
    rectangle('Position',[n- Wbar/2, curve - sem*conf, 2*Wbar/2, sem*conf*2],...
        'EdgeColor',trace,...
        'LineWidth',1);
    hold on
    
    
end

% axes and stuff
ylim([Yinf Ysup]);
set(gca,'FontSize',Font,...
    'XLim',[0 Nbar+1],...
    'XTick',1:Nbar,...    %'YTick',-10:0.2:10,...
    'XTickLabel',varargin);

title(Title);
xlabel(LabelX);
ylabel(LabelY);

end






