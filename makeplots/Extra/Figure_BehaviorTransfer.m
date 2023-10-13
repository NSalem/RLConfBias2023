% plot transfer task behavior for each pair of symbols as well as average
% for each symbol. The behavioral measures are preference, accuracy, 
% confidence and overconfidence

% cd('..')
addpath('ModelingFuncs')
addpath('helperfuncs')

clc;
clear;

%% load data
load('Results\data_all.mat');
addpath('helperfuncs');
whichexp = 1:10;
conf_mat_sess = [];
prevconf_mat_sess = [];
corr_mat_sess = [];
pref = [];
accPost = [];
confpost =[];
ss = [];
prefMat_all= [];
confMat_all = [];
accMat_all = [];
confMatPair_all = [];

pair_mat = [6 5;2 1;8 7;4 3;];
newIdc = pair_mat';

rlVars = load('Results\RLvars_all.mat');
pref  = rlVars.pref;
confpost = rlVars.confpost;
[prefMat_all,accMat_all,confMat_all,confMatPair_all] = makePostMats(rlVars.ss,rlVars.aa,rlVars.confpost);
accPost = squeeze(nanmean(accMat_all./4,3));

pref_col_mat = [
    1,0,0;...
    0,1,0;...
    1,0,1;...
    0,1,1];
col = .65*pref_col_mat([1,1,2,2,3,3,4,4],:);

%%% accuracy and preference all exps %%%%

f1 = figure();

f2 = figure();
for k_meas = 1:2
    switch k_meas
        case 1
            xPost = pref;
            myMat = prefMat_all./4.*100;
            yl = [00,100];
            yLab = 'Preference';
            xLMat = 'S1';
            yLMat = 'S2';
            cL = 'p(choose S1)';
            cLim = [0,100];
        case 2
            xPost = accPost.*100;
            myMat = accMat_all./4.*100;
            yLab = 'Accuracy';
            yl = [0,100];
            xLMat = 'S1';
            yLMat = 'S2';
            cL = 'Accuracy';
            cLim = [0,100];
            
    end
    
    Xx = [xPost(:,pair_mat(1,:)),xPost(:,pair_mat(2,:)),xPost(:,pair_mat(3,:)),xPost(:,pair_mat(4,:))];
    newIdc = pair_mat';
    
    %%% Plot effects per option present %%%
    figure(f1)
    subplot(4,1,k_meas);
    
    dotproperties.useJitter = 1;
    hold on;
    
    plot([0:10],[0:10]*0+50,'--','Color',[.5,.5,.5]);
    
    
    plot([1,2],squeeze(Xx(:,[1,2]))','Color',[0.8,0.8,0.8])
    plot([3,4],squeeze(Xx(:,[3,4]))','Color',[0.8,0.8,0.8])
    plot([5,6],squeeze(Xx(:,[5,6]))','Color',[0.8,0.8,0.8])
    plot([7,8],squeeze(Xx(:,[7,8]))','Color',[0.8,0.8,0.8])
    
    subAvg = mean(Xx,2); %avg per subject accross conditions

    pirateplot(Xx',col,yl(1),yl(2),12,'','',yLab,dotproperties)
    errorbar(8.9,nanmean(subAvg),nanstd(subAvg)./sqrt(numel(subAvg)),'Color','k','Marker','.','MarkerFaceColor','k');
    
    xtickangle(90)
    xlim([0,10]);
    
    figure(f2)
    subplot(4,1,k_meas);
    
    myMatMean = squeeze(nanmean(myMat))';
    myMatMean = myMatMean([newIdc(:)],[newIdc(:)]);
        
    myMatMean = tril(myMatMean);    
    myMatMean(myMatMean==0) = NaN;
    b = imagesc(myMatMean);
    set(b,'AlphaData',~isnan(myMatMean));
    
    lg={'\color{green}G_{75}','\color{green}G_{25}','\color{cyan}G_{75}','\color{cyan}G_{25}',...
        '\color{red}L_{25}','\color{red}L_{75}','\color{magenta}L_{25}','\color{magenta}L_{75}'};
    
    lg={lg{newIdc}};
    
    xtickangle(90);
    xticks(1:8);
    yticks(1:8);
    xticklabels(lg)
    yticklabels(lg);
    xlabel(xLMat,'FontSize',8);ylabel(yLMat,'FontSize',8);
    set(gca,'FontSize',8)
    
    %     xticks([]);yticks([]);
    h = colorbar();
    ylabel(h,cL);
    
    cm = colormap(parula(50));
    caxis(cLim);
    x = linspace(cLim(1),cLim(2),50);
    hold on
end

set(f1,'Position',[0,0,300,900])
set(f2,'Position',[0,0,300,900])

f1 = figure();
f2 = figure();
for k_meas = 1:4
    switch k_meas
        case 1
            xPost = pref(1:90,:,:);
            myMat = prefMat_all(1:90,:,:)./4.*100;
            yl = [00,100];
            yLab = 'Preference';
            xLMat = 'S1';
            yLMat = 'S2';
            cL = 'p(choose S1)';
            cLim = [0,100];
        case 2
            xPost = accPost(1:90,:,:).*100;
            myMat = accMat_all(1:90,:,:)./4.*100;
            yLab = 'Accuracy';
            yl = [0,100];
            xLMat = 'S1';
            yLMat = 'S2';
            cL = 'Accuracy';
            cLim = [0,100];
            
        case 3
            xPost = squeeze(nanmean(confMatPair_all(1:90,:,:),3)).*100;
            myMat = confMat_all(1:90,:,:).*100;
            yl = [50,100];
            yLab = 'Confidence';
            xLMat = 'Unchosen';
            yLMat = 'Chosen';
            cL = 'Confidence';
            cLim = [50,90];
        case 4
            xPost = squeeze(nanmean(confMatPair_all(1:90,:,:),3)).*100-accPost(1:90,:).*100;
            myMat = confMat_all(1:90,:,:).*100-accMat_all(1:90,:,:)./4*100;
            yl = [-20,80];
            yLab = 'Overconfidence';
            xLMat = 'Unchosen';
            yLMat = 'Chosen';
            cL = 'Overconfidence';
            cLim = [-20,80];
            
    end
    
    Xx = [xPost(:,pair_mat(1,:)),xPost(:,pair_mat(2,:)),xPost(:,pair_mat(3,:)),xPost(:,pair_mat(4,:))];
    newIdc = pair_mat';
    
    figure(f1)
    subplot(4,1,k_meas);

    dotproperties.useJitter = 1;
    hold on;
    
    if k_meas<3
        plot([0:10],[0:10]*0+50,'--','Color',[.5,.5,.5]);
    elseif k_meas == 4
        plot([0,10],[0,0],'--','Color',[.5,.5,.5])
    end
    
    plot([1,2],squeeze(Xx(:,[1,2]))','Color',[0.8,0.8,0.8])
    plot([3,4],squeeze(Xx(:,[3,4]))','Color',[0.8,0.8,0.8])
    plot([5,6],squeeze(Xx(:,[5,6]))','Color',[0.8,0.8,0.8])
    plot([7,8],squeeze(Xx(:,[7,8]))','Color',[0.8,0.8,0.8])
    
    subAvg = mean(Xx,2); %avg per subject accross conditions
    pirateplot(Xx',col,yl(1),yl(2),12,'','',yLab,dotproperties)
    errorbar(8.9,nanmean(subAvg),nanstd(subAvg)./sqrt(numel(subAvg)),'Color','k','Marker','.','MarkerFaceColor','k');
    
    xtickangle(90)
    xlim([0,10]);
    
    figure(f2)
    subplot(4,1,k_meas);    
    
    myMatMean = squeeze(nanmean(myMat))';
    myMatMean = myMatMean([newIdc(:)],[newIdc(:)]);
        
    if k_meas <3 %% for preference and accuracy symmetric matrix. For (over)confidence we distinguish chosen/unchosen
        myMatMean = tril(myMatMean);
    end
    
    myMatMean(myMatMean==0) = NaN;
    b = imagesc(myMatMean);
    set(b,'AlphaData',~isnan(myMatMean));
    
    lg={'\color{green}G_{75}','\color{green}G_{25}','\color{cyan}G_{75}','\color{cyan}G_{25}',...
        '\color{red}L_{25}','\color{red}L_{75}','\color{magenta}L_{25}','\color{magenta}L_{75}'};
    
    lg={lg{newIdc}};
    
    xtickangle(90);
    xticks(1:8);
    yticks(1:8);
    xticklabels(lg)
    yticklabels(lg);
    xlabel(xLMat,'FontSize',8);ylabel(yLMat,'FontSize',8);
    set(gca,'FontSize',8)
    
    h = colorbar();
    ylabel(h,cL);
    
    cm = colormap(parula(50));
    caxis(cLim);
    x = linspace(cLim(1),cLim(2),50);
    hold on
end

set(f1,'Position',[0,0,300,900])
set(f2,'Position',[0,0,300,900])


