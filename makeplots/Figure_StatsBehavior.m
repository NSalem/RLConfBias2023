%%% stats and plots of behavior and model predictions across all
%%% experiments

% clc;
clear;
% close all force;

%% load data and model predictions
addpath('helperfuncs')
addpath('helperfuncs/makeLatexTable')

regfile = 'reg_conflogit_learning_dqabs.mat';
% regfile = 'reg_conflogit_norandom_learning_dqabs.mat';
outsuffix = '';

rlvars_all = load('Results\rlvars_all');
% rlvars_all = load('Results\RLvars_norandom_all');
reg = load(['Results\',regfile]);
sims = load('Results\SimulationsSubjParams.mat');


%get variables from data
acc =   squeeze(nanmean(nanmean(rlvars_all.correct,2),4))*100;
conf_1 = squeeze(nanmean(nanmean(reg.confmat,2),4))*100; %conf data

conf = nan(size(acc));
conf(1:90,:) = conf_1;
overconf = conf-acc;

%variables from model fits
imodel = 11;
pcMod =  squeeze(nanmean(nanmean(rlvars_all.pc(imodel,:,:,:,:),3),5))*100;
confBias = erase(reg.formulas,{'conf ~ 1','+ dQabs','+ dQ','+ confprev',' '});
idc = find(reg.isConfPrev & reg.whichLearnModel==imodel & strcmp(confBias,'+Qc'));

confMod_1 = squeeze(nanmean(nanmean(reg.confmodel(:,idc,:,:,:),3),5))*100;
confMod = nan(size(acc));
confMod(1:90,:) = confMod_1;

% variables from simulations
accSims = 100.*squeeze(nanmean(nanmean(sims.outSorted.pc,3),5));
confSims = 100.*squeeze(nanmean(nanmean(nanmean(sims.outConfLSorted(:,:,2,:,:,:),4),6),1));


% make regressors for data and fits
nsub = size(acc,1); %number of participants in data
f1 = [ones(2*(nsub),1);-ones(2*(nsub),1)]; %valence
f2 = repmat([-ones(nsub,1);ones(nsub,1)],2,1); %info
rfx = (1:nsub)';
f3 = repmat(rfx,4,1); %participant
f4 = repmat([rlvars_all.infostr.subexp],1,4)';
% X = [f1,f2,f3,f4];
X = [acc(:),conf(:),confMod(:),f1,f2,double(f3),double(f4)];
% X_conf = X(double(f4)<5,:);
idcconf = find(f3<91);

%make regressors for simulation
nsubs_sim = size(accSims,1);
f1_sims = [ones(2*(nsubs_sim),1);-ones(2*(nsubs_sim),1)]; %valence
f2_sims = repmat([-ones(nsubs_sim,1);ones(nsubs_sim,1)],2,1); %info
rfx_sims = (1:nsubs_sim)';
f3_sims = repmat(rfx_sims,4,1); %participant
%     f4_sims = repmat([rlvars_all.infostr.subexp],1,4)';

% put all exps in one table
varNames = {'Accuracy','Confidence','Overconfidence','AccuracyModel','ConfidenceModel','OverconfidenceModel','Valence','Information','Participant','Experiment'};
tbl = table(acc(:),conf(:),overconf(:),pcMod(:),confMod(:),confMod(:)-pcMod(:),f1(:),f2(:),f3(:),categorical(f4(:)),'VariableNames',varNames);
tbl_confonly = table(acc(idcconf),conf(idcconf),overconf(idcconf),pcMod(idcconf),confMod(idcconf),confMod(idcconf)-pcMod(idcconf),f1(idcconf),f2(idcconf),f3(idcconf),categorical(f4(idcconf)),'VariableNames',varNames);

varNames = {'AccuracySims','ConfidenceSims','OverconfidenceSims','Valence','Information','Participant'};
tbl_sims = table(accSims(:),confSims(:),confSims(:)-accSims(:),f1_sims,f2_sims,f3_sims,'VariableNames',varNames);

lmAll = cell(3,1); %store glm per DV
lmAllPerExp = cell(3,numel(unique(tbl.Experiment))); %store all glms per experiment (DV x Experiment)

lmAllModel = cell(3,1); %store glm per DV
lmAllPerExpModel = cell(3,numel(unique(tbl.Experiment))); %store all glms per experiment (DV x Experiment)


lmAllSims = cell(3,1); %store glm per DV

% loop over dependent variables (acc, conf, overconf) %
measnames = {'Accuracy','Confidence','Overconfidence'};
for imeas = 1:3
    formula = [measnames{imeas},' ~ Valence * Information  + 1 +(1|Participant)'];
    formulaModel = [measnames{imeas},'Model ~ Valence * Information  + 1 +(1|Participant)'];
    formulaSims = [measnames{imeas},'Sims ~ Valence * Information  + 1 + (1|Participant)'];
    formula2 = [measnames{imeas},' ~ Valence * Information  + Experiment+ (1|Participant)'];
    
    if imeas ==1
        lmAllAcc=  fitlme(tbl,formula);
        anovaAllAcc = anova(lmAllAcc,'DFMethod','Satterthwaite');
        lmAllAcc2 = fitlme(tbl,formula2);
        lmAllAccModel=  fitlme(tbl,formulaModel);
        FAllAcc = anova(lmAllAcc2);
    end
    lmAll{imeas} = fitlme(tbl_confonly,formula);
    lmAll2{imeas} = fitlme(tbl_confonly,formula2);
    anovaAll{imeas} = anova(lmAll{imeas},'DFMethod','Satterthwaite');
    lmAllModel{imeas} = fitlme(tbl_confonly,formulaModel);
    lmAllSims{imeas} = fitlme(tbl_sims,formulaSims);
    anovaAllSims{imeas} = anova(lmAllSims{imeas},'DFMethod','Satterthwaite');
    
    for iexp = 1:(numel(unique(tbl.Experiment)))
        if iexp>5 && imeas>1
            continue
        end
        tbl_exp = tbl(double(tbl.Experiment)==iexp,:);
        formula = [measnames{imeas},' ~ Valence * Information + 1  + (1|Participant)'];
        lmAllPerExp{imeas,iexp} = fitlme(tbl_exp,formula);
        formulaModel = [measnames{imeas},'Model ~ Valence * Information + 1  + (1|Participant)'];
        lmAllPerExpModel{imeas,iexp} = fitlme(tbl_exp,formulaModel);
        
    end
end


%% Make matrices of effects over measures, parameters and experiments

%%% change order of experiments: 1-5 no confidence, 6-10 cnfidence
expid = [6:10,1:5];

%%% plot coefficients
parsExp = nan(3,4,12);
parsExpSE = nan(3,4,12);
parsExpModel = nan(3,4,12);
parsExpModelSE = nan(3,4,12);
parsExpSims = nan(3,4,12);
parsExpSimsSE = nan(3,4,12);

tExp = nan(3,4,12);
tExpModel = nan(3,4,12);
tExpSims = nan(3,4,12);

for imeas = 1:3
    for ipar = 1:4
        for iexp = 1:12
            if (iexp ==11)
                
                [b,~,stats] = fixedEffects(lmAll{imeas},'DFMethod','Satterthwaite');
                parsExp(imeas,ipar,iexp) = stats.Estimate(ipar);
                parsExpSE(imeas,ipar,iexp) = stats.SE(ipar);
                tExp(imeas,ipar,iexp) = stats.tStat(ipar);
                
                [b,~,stats] = fixedEffects(lmAllModel{imeas},'DFMethod','Satterthwaite');
                parsExpModel(imeas,ipar,iexp) = stats.Estimate(ipar);
                parsExpModelSE(imeas,ipar,iexp) = stats.SE(ipar);
                tExpModel(imeas,ipar,iexp) = stats.tStat(ipar);
                
                [b,~,stats] = fixedEffects(lmAllSims{imeas},'DFMethod','Satterthwaite');
                parsExpSims(imeas,ipar,iexp) = stats.Estimate(ipar);
                parsExpSimsSE(imeas,ipar,iexp) = stats.SE(ipar);
                tExpSims(imeas,ipar,iexp) =stats.tStat(ipar);
                
            elseif (iexp ==12 && imeas==1)
                
                [b,~,stats] = fixedEffects(lmAllAcc,'DFMethod','Satterthwaite');
                parsExp(imeas,ipar,iexp) = stats.Estimate(ipar);
                parsExpSE(imeas,ipar,iexp) = stats.SE(ipar);
                tExp(imeas,ipar,iexp) = stats.tStat(ipar);
                
                [b,~,stats] = fixedEffects(lmAllModel{imeas},'DFMethod','Satterthwaite');
                parsExpModel(imeas,ipar,iexp) = stats.Estimate(ipar);
                parsExpModelSE(imeas,ipar,iexp) = stats.SE(ipar);
                tExpModel(imeas,ipar,iexp) = stats.tStat(ipar);
                
            elseif ~(iexp>5 && imeas>1)
                thisexp = expid(iexp);
                [b,~,stats] = fixedEffects(lmAllPerExp{imeas,iexp},'DFMethod','Satterthwaite');
                parsExp(imeas,ipar,thisexp) = stats.Estimate(ipar);
                parsExpSE(imeas,ipar,thisexp) = stats.SE(ipar);
                tExp(imeas,ipar,thisexp) = stats.tStat(ipar);
                
                [b,~,stats] = fixedEffects(lmAllPerExpModel{imeas,iexp},'DFMethod','Satterthwaite');
                parsExpModel(imeas,ipar,thisexp) = stats.Estimate(ipar);
                parsExpModelSE(imeas,ipar,thisexp) = stats.SE(ipar);
                tExpModel(imeas,ipar,thisexp) = stats.tStat(ipar);
                
            end
        end
    end
end

%%  Plot effects
coldat = [.2,.2,.2];
colmod = [0.4,0.4,0.9];
colsims = [0,0.7,0.1];
figure()
parnames = {'B_0','B_V','B_I','B_{V:I}'};
for imeas = 1:3
    for ipar = 1:4
        subplot(3,4,sub2ind([4,3],ipar,imeas))
        hold on
        
        ylabel(parnames{ipar})
        l1 =errorbar([1:12]-0.1,squeeze(parsExp(imeas,ipar,:)),squeeze(parsExpSE(imeas,ipar,:)),'Marker','None','LineStyle','None','LineWidth',1.5,'Color',coldat);
        
        title(['Learning: ',measnames{imeas}]);
        
        l2 = errorbar([1:12],squeeze(parsExpModel(imeas,ipar,:)),...
            squeeze(parsExpModelSE(imeas,ipar,:)),'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colmod);
        
        l3 = errorbar([1:12]+0.1,squeeze(parsExpSims(imeas,ipar,:)),...
            squeeze(parsExpSimsSE(imeas,ipar,:)),'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colsims);
        
        if numel(unique(sign(get(gca,'YLim'))))>1
            hl = plot([0,13],[0,0],'k--')
            uistack(hl,'bottom')
        end
        
        yl = get(gca,'YLim');
        
        mymax = max([squeeze(parsExp(imeas,ipar,:))+squeeze(parsExpSE(imeas,ipar,:));squeeze(parsExpModel(imeas,ipar,:))+squeeze(parsExpModelSE(imeas,ipar,:))]);
        mymin = min([squeeze(parsExp(imeas,ipar,:))+squeeze(parsExpSE(imeas,ipar,:));squeeze(parsExpModel(imeas,ipar,:))+squeeze(parsExpModelSE(imeas,ipar,:))]);
        
        for iexp = 1:12
            if iexp >= 11
                w = 'bold';
            else
                w = 'normal';
            end
            
            if abs(tExp(imeas,ipar,iexp))>1.96
                text(iexp,mymax+0.2*(yl(2)-yl(1)),'*','Color',coldat,'FontWeight',w)
            end
            if abs(tExpModel(imeas,ipar,iexp))>1.96
                text(iexp,mymax+0.1*(yl(2)-yl(1)),'*','Color',colmod,'FontWeight',w)
            end
        end
        
        if imeas ==1 && (ipar) ==1
            ylim([mymin-0.8*diff(yl),mymax+0.25*(diff(yl))])
        else
            ylim([yl(1),mymax+0.25*(diff(yl))])
        end
        
        yl = get(gca,'YLim');
        sq = fill([10.5,13,13,10.5],[yl(1),yl(1),yl(2),yl(2)],[0.5,0.5,0.5],'FaceAlpha',0.25,'EdgeAlpha',0)
        uistack(sq,'bottom')
        
        vl = plot([10.5,10.5],yl,'k-')
        
        xticks(1:12)
        xticklabels({'1','2','3','4','5','6','7','8','9','10','Conf','All'});
        %             xlabel('Experiment');
        xtickangle(90)
        xlim([0,13])
        if imeas ==1 && ipar ==1
            leg = legend([l1,l2,l3],{'Data','Predicted','Simulation'},'Location','best');
            leg.BoxFace.ColorType='truecoloralpha';
            leg.BoxFace.ColorData=(uint8(255*[1 1 1 0.75]'));
        end
        
        set(gca,'FontName','Arial','FontSize',11)
    end
end
set(gcf,'Position',[300,100,1200,700])
saveas(gcf,['Plots/coeffLT',outsuffix,'.png']);
saveas(gcf,['Plots/coeffLT',outsuffix,'.svg']);


%%  Plot sorted coefficients (all exp together)

coldat = [.2,.2,.2];
colmod = [0.4,0.4,0.9];
colsims = [0,0.7,0.1];
figure()
plotMatM = zeros(5,numel(measnames),4);
plotMatSE = zeros(5,numel(measnames),4);
plotMatM(1,:,:) = parsExp(:,:,12);
plotMatSE(1,:,:) = parsExpSE(:,:,12);
plotMatM(2,:,:) = parsExpModel(:,:,12);
plotMatSE(2,:,:) = parsExpModelSE(:,:,12);
parnames = {'0','V','I','V:I'};
% parnames = {'\beta_{0}','\beta_{Valence}','\beta_{Information}','\beta_{Interaction}'};
for imeas = 1
    subplot(4,1,imeas)
    hold on
    %sort by effect size
    dumM = parsExp(imeas,2:4,12);
    [~,id] = sort((dumM),'descend');
    dumM = dumM(id);
    dumSE = parsExpSE(imeas,2:4,12);
    dumSE = dumSE(id);
    
    offset = 0;
    l1 = errorbar([1:3]-0.1,dumM,squeeze(dumSE),'Marker','None','LineStyle','None','LineWidth',1.5,'Color',coldat);
    title(['Learning: ', measnames{imeas},'(all exp.)']);
    ylabel('Regression fixed-effect')
    yl = get(gca,'YLim');
    mymax = max(max(squeeze(plotMatM(:,imeas,2:4))+squeeze(plotMatSE(:,imeas,2:4))));
    mymin = min(min(squeeze(plotMatM(:,imeas,2:4))-squeeze(plotMatSE(:,imeas,2:4))));
    
    dummodelM = parsExpModel(imeas,2:4,12);
    dummodelSE = parsExpModelSE(imeas,2:4,12);
    dummodelM = squeeze(dummodelM(id));
    dummodelSE = squeeze(dummodelSE(id));
    
    l2 = errorbar([1:3],dummodelM,dummodelSE,'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colmod);
    
    dumsimsM = parsExpSims(imeas,2:4,11);
    dumsimsSE = parsExpSimsSE(imeas,2:4,11);
    dumsimsM = squeeze(dumsimsM(id));
    dumsimsSE = squeeze(dumsimsSE(id));
    
    l3 = errorbar([1:3]+0.1,dumsimsM,dumsimsSE,'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colsims);
    
    %show significance
    dumtExp = tExp(imeas,2:end,12);
    dumtExp = dumtExp(id);
    dumtExpModel = squeeze(tExpModel(imeas,2:end,12));
    dumtExpModel = dumtExpModel(id);
    
    for ipar = 1:numel(id)
        if abs(dumtExp(ipar))>1.96
            text(ipar-offset,mymax+0.2*(yl(2)-yl(1)),'*','Color',coldat,'HorizontalAlignment','center')
        end
        
        if abs(dumtExpModel(ipar))>1.96
            text(ipar+offset,mymax+0.1*(yl(2)-yl(1)),'*','Color',colmod,'HorizontalAlignment','center')
        end
    end
    
    yl = get(gca,'YLim');
    
    if imeas ==1 && (ipar) ==1
        ylim([mymin-0.8*diff(yl),mymax+0.25*(diff(yl))])
    else
        ylim([yl(1),mymax+0.25*(diff(yl))])
    end
    %
    yl = get(gca,'YLim');
    if numel(unique(sign(get(gca,'YLim'))))>1
        plot([0,4],[0,0],'k--')
        %         uistack(hl,'bottom')
    end
    
    xticks(1:3)
    dumxlabels = parnames(2:4);
    xticklabels(dumxlabels(id))
    
    xtickangle(90)
    xlim([0,numel(dumxlabels)+1])
    if imeas ==1
        leg = legend([l1,l2,l3],{'Data','Predicted','Simulation'},'Location','south');
        leg.BoxFace.ColorType='truecoloralpha';
        leg.BoxFace.ColorData=(uint8(255*[1 1 1 0.75]'));
    end
    set(gca,'FontName','Arial','FontSize',11)
end

set(gcf,'Position',[0,0,200,1000])
saveas(gcf,'Plots/coefSortedLT_allData.svg');

%%% Only confidence exps
figure()
plotMatM = zeros(5,numel(measnames),4);
plotMatSE = zeros(5,numel(measnames),4);
plotMatM(1,:,:) = parsExp(:,:,11);
plotMatSE(1,:,:) = parsExpSE(:,:,11);
plotMatM(2,:,:) = parsExpModel(:,:,11);
plotMatSE(2,:,:) = parsExpModelSE(:,:,11);
parnames = {'0','V','I','V:I'};
% parnames = {'\beta_{0}','\beta_{Valence}','\beta_{Information}','\beta_{Interaction}'};
for imeas = 1:numel(measnames)
    subplot(4,1,imeas)
    hold on
    %sort by effect size
    dumM = parsExp(imeas,2:4,11);
    [~,id] = sort((dumM),'descend');
    dumM = dumM(id);
    dumSE = parsExpSE(imeas,2:4,11);
    dumSE = dumSE(id);
    
    offset = 0;
    l1 = errorbar([1:3]-0.1,dumM,squeeze(dumSE),'Marker','None','LineStyle','None','LineWidth',1.5,'Color',coldat);
    title(['Learning: ', measnames{imeas}]);
    ylabel('Regression fixed-effect')
    yl = get(gca,'YLim');
    mymax = max(max(squeeze(plotMatM(:,imeas,2:4))+squeeze(plotMatSE(:,imeas,2:4))));
    mymin = min(min(squeeze(plotMatM(:,imeas,2:4))-squeeze(plotMatSE(:,imeas,2:4))));
    
    dummodelM = parsExpModel(imeas,2:4,11);
    dummodelSE = parsExpModelSE(imeas,2:4,11);
    dummodelM = squeeze(dummodelM(id));
    dummodelSE = squeeze(dummodelSE(id));
    
    l2 = errorbar([1:3],dummodelM,dummodelSE,'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colmod);
    
    dumsimsM = parsExpSims(imeas,2:4,11);
    dumsimsSE = parsExpSimsSE(imeas,2:4,11);
    dumsimsM = squeeze(dumsimsM(id));
    dumsimsSE = squeeze(dumsimsSE(id));
    
    l3 = errorbar([1:3]+0.1,dumsimsM,dumsimsSE,'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colsims);
    
    %show significance
    dumtExp = tExp(imeas,2:end,11);
    dumtExp = dumtExp(id);
    dumtExpModel = squeeze(tExpModel(imeas,2:end,11));
    dumtExpModel = dumtExpModel(id);
    
    for ipar = 1:numel(id)
        if abs(dumtExp(ipar))>1.96
            text(ipar-offset,mymax+0.2*(yl(2)-yl(1)),'*','Color',coldat,'HorizontalAlignment','center')
        end
        
        if abs(dumtExpModel(ipar))>1.96
            text(ipar+offset,mymax+0.1*(yl(2)-yl(1)),'*','Color',colmod,'HorizontalAlignment','center')
        end
    end
    
    yl = get(gca,'YLim');
    
    if imeas ==1 && (ipar) ==1
        ylim([mymin-0.8*diff(yl),mymax+0.25*(diff(yl))])
    else
        ylim([yl(1),mymax+0.25*(diff(yl))])
    end
    %
    yl = get(gca,'YLim');
    if numel(unique(sign(get(gca,'YLim'))))>1
        plot([0,4],[0,0],'k--')
    end
    
    xticks(1:3)
    dumxlabels = parnames(2:4);
    xticklabels(dumxlabels(id))
    
    xtickangle(90)
    xlim([0,numel(dumxlabels)+1])
    if imeas ==1
        leg = legend([l1,l2,l3],{'Data','Predicted','Simulation'},'Location','best');
        leg.BoxFace.ColorType='truecoloralpha';
        leg.BoxFace.ColorData=(uint8(255*[1 1 1 0.75]'));
    end
    set(gca,'FontName','Arial','FontSize',11)
end

set(gcf,'Position',[0,0,200,1000])
saveas(gcf,'Plots/coefSortedLT_confData.svg');

%% Write latex output
%%% -iterate over measures
%%% - save single tex file with all tables
latexCode = [];
for imeas = 1:3
    if imeas ==1
        %%%% make table with all data for accuracy
        thisTable =  convert_LMM2latex(...
            {lmAllAcc,lmAllAccModel},... % array of glme results
            {'Data','Model'},... % name of the models
            {'tStat'},... % information required
            2,... % round digits
            ['reg',measnames{imeas},'_allexp'],... % table label
            'vertical'); % orientation (vertical or rotated)
        
        thisTable=replaceWords(thisTable,...
            {'Insert Title Here.'},... % old names
            {['Regression results for ',measnames{imeas}, ' in Learning Task (all data)']}); % % new names
        
        latexCode =[latexCode;thisTable];
        
    end
    
    thisTable =  convert_LMM2latex(...
        {lmAll{imeas},lmAllModel{imeas},lmAllSims{imeas}},... % array of glme results
        {'Data','Model','Simulations'},... % name of the models
        {'tStat'},... % information required
        2,... % round digits
        ['reg',measnames{imeas},'_allexp'],... % table label
        'vertical'); % orientation (vertical or rotated)
    
    thisTable=replaceWords(thisTable,...
        {'Insert Title Here.'},... % old names
        {['Regression results for ',measnames{imeas}, ' in Learning Task (confidence experiments)']}); % % new names
    
    latexCode =[latexCode;thisTable];
    
    thisTable = convert_LMM2latex(...
        lmAllPerExp(imeas,1:5),... % array of glme results
        {'Exp. 1', 'Exp. 2', 'Exp. 3', 'Exp. 4','Exp.5'},... % name of the models
        {'tStat'},... % information required
        2,... % round digits
        ['reg',measnames{imeas},'_confexps'],... % table label
        'vertical'); % orientation (vertical or rotated)
    
    thisTable=replaceWords(thisTable,...
        {'Insert Title Here.'},... % old names
        {['Regression results for ',measnames{imeas}, ' for each confidence experiment']}); % % new names
    
    latexCode =[latexCode;thisTable];
    
    if imeas ==1
        thisTable = convert_LMM2latex(...
            lmAllPerExp(imeas,6:10),... % array of glme results
            {'Exp. 6', 'Exp. 7', 'Exp. 8', 'Exp. 9','Exp.10'},... % name of the models
            {'tStat'},... % information required
            2,... % round digits
            ['reg',measnames{imeas},'_noconfexps'],... % table label
            'vertical'); % orientation (vertical or rotated)
        
        thisTable=replaceWords(thisTable,...
            {'Insert Title Here.'},... % old names
            {['Regression results for ',measnames{imeas}, ' for each non-confidence experiment']}); % % new names
        
        latexCode =[latexCode;thisTable];
        
    end
    
    
    thisTable = convert_LMM2latex(...
        lmAllPerExpModel(imeas,1:5),... % array of glme results
        {'Exp. 1', 'Exp. 2', 'Exp. 3', 'Exp. 4','Exp.5'},... % name of the models
        {'tStat'},... % information required
        2,... % round digits
        ['reg',measnames{imeas},'Model_confexps'],... % table label
        'vertical'); % orientation (vertical or rotated)
    
    thisTable=replaceWords(thisTable,...
        {'Insert Title Here.'},... % old names
        {['Regression results for model-predicted ',measnames{imeas}, ' for each confidence experiment']}); % % new names
    
    latexCode =[latexCode;thisTable];
    
    
    if imeas ==1
        thisTable = convert_LMM2latex(...
            lmAllPerExpModel(imeas,6:10),... % array of glme results
            {'Exp. 6', 'Exp. 8', 'Exp. 8', 'Exp. 9','Exp.10'},... % name of the models
            {'tStat'},... % information required
            2,... % round digits
            ['reg',measnames{imeas},'Model_noconfexps'],... % table label
            'vertical'); % orientation (vertical or rotated)
        
        thisTable=replaceWords(thisTable,...
            {'Insert Title Here.'},... % old names
            {['Regression results for model-predicted ',measnames{imeas}, ' for each non-confidence experiment']}); % % new names
        
        latexCode =[latexCode;thisTable];
        
    end
end

%%% write file
fid=fopen('Results/LatexTables/statsLT.tex','w'); fprintf(fid,'%s\n',latexCode{:}); fclose(fid);

