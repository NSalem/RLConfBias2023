% Plot behavior and model fits (based
% on best learning+confidence model) for Transfer Task. 
% Plots incluided: Violin plots and option pair matrices for preference, 
% accuracy, confidence, overconfidence.

% cd('..')
addpath('ModelingFuncs')
addpath('helperfuncs')

clc;
clear;
% close all force;

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

reg = load('Results\reg_conflogit_posttest_dqabs.mat');
rlVars = load('Results\RLvars_all.mat');

ilearnmodel = [11];

iconfmodel = find(reg.whichLearnModel == ilearnmodel & strcmp(reg.confBias, '+Qc') & ~reg.isConfPrev);

pref  = rlVars.pref;
confpost = rlVars.confpost;
[prefMat_all,accMat_all,confMat_all,confMatPair_all] = makePostMats(rlVars.ss,rlVars.aa,rlVars.confpost);
accPost = squeeze(nanmean(accMat_all./4,3));


%% simulate choices model
accModMat = NaN(size(rlVars.Q,2),8,8);
prefModMat =  NaN(size(rlVars.Q,2),8,8);
QfAll = NaN(size(rlVars.Q,2),8);
for isub = 1:size(rlVars.Q,2)
    %%% create Q-learner w/subj beta (don't need learning params for choice
    %%% only
    paramstruct = struct('beta',rlVars.params_allexp{isub,ilearnmodel}(1));
    l = qLearner(paramstruct);
    
    %%get final q-values, reorder to PL75,PL25,PG25,PG75,CL75,CL25,CG25,CG75
    nsess = rlVars.infostr.sub_nsess(isub);
    ntrl = rlVars.infostr.sub_ntrials(isub);
    Qf = squeeze(rlVars.Q(ilearnmodel,isub,nsess,:,:,ntrl));
    id = [3,7,1,5,4,8,2,6];
    Qf = Qf(id);
    QfAll(isub,:) = Qf;
    %%% get real utility
    optP = [0.75,0.25,0.75,0.25,0.25,0.75,0.25,0.75];
    optU = [optP*1+(1-optP)*0.1]; %option utility (values in exp 1 are different but same order)
    optU(5:8) = -optU(5:8);
    optU = optU([8,7,2,1,6,5,4,3]);
    %%% for each option(except L75), compute model expected accuracy (choice rate vs worse options)
    p = nan(8);
    for iopt = 1:8
        for iopt2 = 1:8
            optDiff = optU(iopt)-optU(iopt2);
            [a,p] = l.chooseActionPostLearning(Qf(iopt2),Qf(iopt));
            if a ~= 2
                p = 1-p; %p is prob of a, make it prob of option 1
            end
            prefModMat(isub,iopt,iopt2) = p;
            if optDiff<0 %if option 2 better, p(correct) = p(2)
                p = 1-p;
            end
            if optDiff~=0
                accModMat(isub,iopt,iopt2) = p;
            end
        end
    end
end

accMod = squeeze(nanmean(accModMat,3));
prefMod = squeeze(rlVars.ppc(ilearnmodel,:,:));
prefMod = prefMod(:,pair_mat');

pref_col_mat = [
    1,0,0;...
    0,1,0;...
    1,0,1;...
    0,1,1];
col = .65*pref_col_mat([1,1,2,2,3,3,4,4],:);


%% accuracy and preference all exps 
f1 = figure();
f2 = figure();
for k_meas = 1:2
    switch k_meas
        case 1
            xPost = pref;
            xMod = prefMod.*100;
            myMat = prefMat_all./4.*100;
            modMat = prefModMat.*100;
            yl = [00,100];
            yLab = 'Preference';
            xLMat = 'S1';
            yLMat = 'S2';
            cL = 'p(choose S1)';
            cLim = [0,100];
        case 2
            xPost = accPost.*100;
            xMod = accMod.*100;
            myMat = accMat_all./4.*100;
            modMat = accModMat.*100;
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
    %     subplot(4,2,sub2ind([2,4],1,k_meas));
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
    %         pirateplot(Xx(:,[2:4,6:8])',col([2:4,6:8],:),yl(1),yl(2),12,'','Chosen option',yLab,dotproperties)
    %         pirateplot(Xx(:,6:8)',col(5:8,:),0,100,12,'Data','Chosen option','Preference (%)',dotproperties)
    pirateplot(Xx',col,yl(1),yl(2),12,'','',yLab,dotproperties)
    %     xticklabels([fliplr({'PG_7_5','PG_2_5','PL_2_5','PL_7_5'}),fliplr({'CG_7_5','CG_2_5','CL_2_5','CL_7_5'})])
    %         xticklabels([fliplr({'PG_7_5','PG_2_5','PL_2_5'}),fliplr({'CG_7_5','CG_2_5','CL_2_5'})])
    errorbar(8.9,nanmean(subAvg),nanstd(subAvg)./sqrt(numel(subAvg)),'Color','k','Marker','.','MarkerFaceColor','k');
    
    for iopt = 1:8
        m = nanmean(xMod(:,iopt));
        se = nanstd(xMod(:,iopt))./sqrt(size(xMod,1));
        errorbar(iopt,m,se,'MarkerFaceColor','white','Color','k','Marker','o')
    end
    modelSubAvg = nanmean(xMod,2);
    mAll = nanmean(modelSubAvg);
    seAll = nanstd(modelSubAvg)./sqrt(size(modelSubAvg,1));
    
    errorbar(9.1,mAll,seAll,'Color','k','Marker','d','MarkerFaceColor','white');
    
    xtickangle(90)
    xlim([0,10]);
    
    %     subplot(4,2,sub2ind([2,4],2,k_meas));
    figure(f2)
    subplot(4,1,k_meas);
    
    myMatMean = squeeze(nanmean(myMat))';
    myMatMean = myMatMean([newIdc(:)],[newIdc(:)]);
    
    
    modelMatMean = squeeze(nanmean(modMat))';
    
    myMatMean = tril(myMatMean);
    modelMatMean = tril(modelMatMean);
    
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
    % plot(1,2,confMatMean(1,2))
    caxis(cLim);
    x = linspace(cLim(1),cLim(2),50);
    hold on
    for ir = 1:8
        for ic = 1:8
            if ~isnan(myMatMean(ir,ic))
                [~,id]=min(abs(modelMatMean(ir,ic)-x));
                scatter(ic,ir,50,'MarkerFaceColor',cm(id,:),'Marker','o','MarkerEdgeColor',[0,0,0]);
            end
        end
    end
end

set(f1,'Position',[0,0,300,900])
set(f2,'Position',[0,0,300,900])

% saveas(f1,'Plots/modelFitsTransferAllData_Matrices.svg')
% saveas(f2,'Plots/modelFitsTransferAllData_Violins.svg')

%% conf exps only %%
confMod = squeeze(reg.confmodelpostPresent(:,iconfmodel,:)).*100;
confMod = confMod(:,pair_mat');
confModMat = squeeze(reg.confmodelmatChoice(iconfmodel,:,newIdc,newIdc)).*100;

f1 = figure();
f2 = figure();
for k_meas = 1:4
    switch k_meas
        case 1
            xPost = pref(1:90,:,:);
            xMod = prefMod(1:90,:,:).*100;
            myMat = prefMat_all(1:90,:,:)./4.*100;
            modMat = prefModMat(1:90,:,:).*100;
            yl = [00,100];
            yLab = 'Preference';
            xLMat = 'S1';
            yLMat = 'S2';
            cL = 'p(choose S1)';
            cLim = [0,100];
        case 2
            xPost = accPost(1:90,:,:).*100;
            xMod = accMod(1:90,:,:).*100;
            myMat = accMat_all(1:90,:,:)./4.*100;
            modMat = accModMat(1:90,:,:).*100;
            yLab = 'Accuracy';
            yl = [0,100];
            xLMat = 'S1';
            yLMat = 'S2';
            cL = 'Accuracy';
            cLim = [0,100];
            
        case 3
            xPost = squeeze(nanmean(confMatPair_all(1:90,:,:),3)).*100;
            xMod = confMod(1:90,:);
            myMat = confMat_all(1:90,:,:).*100;
            modMat = confModMat(1:90,:,:);
            yl = [50,100];
            yLab = 'Confidence';
            xLMat = 'Unchosen';
            yLMat = 'Chosen';
            cL = 'Confidence';
            cLim = [50,90];
        case 4
            xPost = squeeze(nanmean(confMatPair_all(1:90,:,:),3)).*100-accPost(1:90,:).*100;
            xMod = confMod - accMod(1:90,:).*100;
            myMat = confMat_all(1:90,:,:).*100-accMat_all(1:90,:,:)./4*100;
            modMat = confModMat(1:90,:,:) - accModMat(1:90,:,:).*100;
            yl = [-20,80];
            yLab = 'Overconfidence';
            xLMat = 'Unchosen';
            yLMat = 'Chosen';
            cL = 'Overconfidence';
            cLim = [-20,80];
            
    end
    
    Xx = [xPost(:,pair_mat(1,:)),xPost(:,pair_mat(2,:)),xPost(:,pair_mat(3,:)),xPost(:,pair_mat(4,:))];
    newIdc = pair_mat';
    
    %%% Plot effects per option present %%%
%     subplot(4,2,sub2ind([2,4],1,k_meas));
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
    %         pirateplot(Xx(:,[2:4,6:8])',col([2:4,6:8],:),yl(1),yl(2),12,'','Chosen option',yLab,dotproperties)
    %         pirateplot(Xx(:,6:8)',col(5:8,:),0,100,12,'Data','Chosen option','Preference (%)',dotproperties)
    pirateplot(Xx',col,yl(1),yl(2),12,'','',yLab,dotproperties)
    %     xticklabels([fliplr({'PG_7_5','PG_2_5','PL_2_5','PL_7_5'}),fliplr({'CG_7_5','CG_2_5','CL_2_5','CL_7_5'})])
    %         xticklabels([fliplr({'PG_7_5','PG_2_5','PL_2_5'}),fliplr({'CG_7_5','CG_2_5','CL_2_5'})])
    errorbar(8.9,nanmean(subAvg),nanstd(subAvg)./sqrt(numel(subAvg)),'Color','k','Marker','.','MarkerFaceColor','k');
    
    for iopt = 1:8
        m = nanmean(xMod(:,iopt));
        se = nanstd(xMod(:,iopt))./sqrt(size(xMod,1));
        errorbar(iopt,m,se,'MarkerFaceColor','white','Color','k','Marker','o')
    end
    modelSubAvg = nanmean(xMod,2);
    mAll = nanmean(modelSubAvg);
    seAll = nanstd(modelSubAvg)./sqrt(size(modelSubAvg,1));
    
    errorbar(9.1,mAll,seAll,'Color','k','Marker','d','MarkerFaceColor','white');
    
    xtickangle(90)
    xlim([0,10]);
    
    figure(f2)
    subplot(4,1,k_meas);    
    
    myMatMean = squeeze(nanmean(myMat))';
    myMatMean = myMatMean([newIdc(:)],[newIdc(:)]);
    
    
    modelMatMean = squeeze(nanmean(modMat))';
    
    if k_meas <3 %% for preference and accuracy symmetric matrix. For (over)confidence we distinguish chosen/unchosen
        myMatMean = tril(myMatMean);
        modelMatMean = tril(modelMatMean);
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
    
    %     xticks([]);yticks([]);
    h = colorbar();
    ylabel(h,cL);
    
    cm = colormap(parula(50));
    % plot(1,2,confMatMean(1,2))
    caxis(cLim);
    x = linspace(cLim(1),cLim(2),50);
    hold on
    for ir = 1:8
        for ic = 1:8
            if ~isnan(myMatMean(ir,ic))
                [~,id]=min(abs(modelMatMean(ir,ic)-x));
                scatter(ic,ir,50,'MarkerFaceColor',cm(id,:),'Marker','o','MarkerEdgeColor',[0,0,0]);
            end
        end
    end
end

set(f1,'Position',[0,0,300,900])
set(f2,'Position',[0,0,300,900])

% saveas(f1,'Plots/modelFitsTransferConfData_Matrices.svg')
% saveas(f1,'Plots/modelFitsTransferConfData_Violins.svg')


