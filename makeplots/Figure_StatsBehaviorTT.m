%% Transfer task

%% load data %
% load('Results\data_all.mat');
addpath('helperfuncs');
addpath('ModelingFuncs\')
addpath('helperfuncs/makeLatexTable')

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

regfile = 'reg_conflogit_posttest_dqabs.mat';
outsuffix = '';

rlvars_all = load('Results\rlvars_all');
reg = load(['Results\',regfile]);
sims = load('Results\SimulationsSubjParams.mat');

ilearnmodel = [11];
confBias = erase(reg.formulas,{'conf ~ 1','+ dQabs','+ dQ','+ confprev',' '});
iconfmodel = find(~reg.isConfPrev & reg.whichLearnModel==ilearnmodel & strcmp(confBias,'+Qc'));


%% simulate choices model
pref  = rlvars_all.pref;
confpost = rlvars_all.confpost;

%%% get matrices for data
[prefMat_all,accMat_all,confMat_all,confMatPair_all] = makePostMats(rlvars_all.ss,rlvars_all.aa,rlvars_all.confpost);

accPost = squeeze(nanmean(accMat_all./4,3));
accModMat = NaN(size(rlvars_all.Q,2),8,8);
prefModMat =   NaN(size(rlvars_all.Q,2),8,8);
QfAll = NaN(size(rlvars_all.Q,2),8);

for isub = 1:size(rlvars_all.Q,2)
    %%% create Q-learner w/subj beta (don't need learning params for choice
    %%% only
    paramstruct = struct('beta',rlvars_all.params_allexp{isub,ilearnmodel}(find(strcmp(rlvars_all.modelsinfo{ilearnmodel}.paramnames,'beta'))));
    paramstruct.betaTT = rlvars_all.params_allexp{isub,ilearnmodel}(find(strcmp(rlvars_all.modelsinfo{ilearnmodel}.paramnames,'betaTT')));

    l = qLearner(paramstruct);
    %%get final q-values, reorder to PL75,PL25,PG25,PG75,CL75,CL25,CG25,CG75
    nsess = rlvars_all.infostr.sub_nsess(isub);
    ntrl = rlvars_all.infostr.sub_ntrials(isub);
    Qf = squeeze(rlvars_all.Q(ilearnmodel,isub,nsess,:,:,ntrl));
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
    if all(isnan(squeeze(accModMat(isub,:))))
        keyboard()
    end
end

accMod = squeeze(nanmean(accModMat,3));
prefMod = [squeeze(rlvars_all.ppc(ilearnmodel,:,:))];
prefMod = prefMod(:,pair_mat');
confMod_1 = squeeze(reg.confmodelpostPresent(:,iconfmodel,:));
confMod_1 = confMod_1(:,pair_mat').*100;

confModMat_1 = squeeze(reg.confmodelmatChoice(iconfmodel,:,newIdc,newIdc))*100;

confMod = nan(size(accMod));
confMod(1:90,:) = confMod_1;

confModMat = nan(size(accModMat));
confModMat(1:90,:,:) = confModMat_1;

accSims = 100.*squeeze(nanmean(sims.accPost));
confSims = 100.*squeeze(nanmean(sims.confPost(:,2,:,:)));

% prefSims = squeeze(nanmean(sims.pref));
prefSimsMat =squeeze(nanmean(sims.prefMat))./4*100;
prefSims = nanmean(prefSimsMat,3);
accSimsMat = squeeze(nanmean(sims.accMat))./4.*100;
confSimsMat= 100.*squeeze(nanmean(sims.confChoicePairMat(:,2,:,:,:)))

%% do stats over each DV
figure()
set(gcf,'Position',[0,0,600,900])
lmAll = cell(4,1);
lmAllExp = cell(4,10);

lmAllModel = cell(4,1);
lmAllExpModel = cell(4,10);
latexCode = [];
for k_meas = 1:4
    switch k_meas
        case 1
            xPost = pref;
            xMod = prefMod.*100;
            xSims = prefSims;
            myMat = prefMat_all./4.*100;
            modMat = prefModMat.*100;
            simsMat = prefSimsMat.*100;
        case 2
            xPost = accPost.*100;
            xMod = accMod.*100;
            xSims = accSims;
            myMat = accMat_all./4.*100;
            modMat = accModMat.*100;
            simsMat = accSimsMat;
        case 3
            xPost = squeeze(nanmean(confMatPair_all(:,:,:),3)).*100;
            xMod = confMod;
            xSims = confSims;
            myMat = confMat_all(:,:,:).*100;
            modMat = confModMat(:,:,:);
            simsMat = confSimsMat;
 
        case 4
            xPost =squeeze(nanmean(confMatPair_all(:,:,:),3)).*100-accPost(:,:).*100;
            xMod = confMod(:,:)- accMod(:,:).*100;
            xSims = confSims-accSims;
            myMat = confMat_all(:,:,:).*100-accMat_all(:,:,:)./4*100;
            modMat = confModMat(:,:,:) - accModMat(:,:,:).*100;
            simsMat = confSimsMat - accSimsMat;
    end
    
    Xx = [xPost(:,pair_mat(1,:)),xPost(:,pair_mat(2,:)),xPost(:,pair_mat(3,:)),xPost(:,pair_mat(4,:))];
    XxMod = xMod;
    XxSims = xSims(:,newIdc);
    
    %%% Make regressors for data and fits %%%
    f1 = repmat([-1,-1,1,1,-1,-1,1,1],size(XxMod,1),1);
    f1 = f1(:);
    f2 = repmat([-1,1,-1,1,-1,1,-1,1],size(XxMod,1),1); %correctness in original context
    f2 = f2(:);
    f3 = repmat([-ones(1,4),ones(1,4)],size(XxMod,1),1);
    f3 = f3(:);
    rfx = (1:size(XxMod,1))';
    f4 = repmat(rfx,8,1);
    f5 = repmat([rlvars_all.infostr.subexp],1,8)';
    
    X = [f1,f2,f3,f4,f5];
    vnames = {'Y','Valence','Correctness','Information','Participant','Experiment'};
    vnamesSims = {'Y','Valence','Correctness','Information','Participant'};
    
    %%% Make regressors for simulations %%%
    f1_sims = repmat([-1,-1,1,1,-1,-1,1,1],size(XxSims,1),1);
    f1_sims = f1_sims(:);
    f2_sims = repmat([-1,1,-1,1,-1,1,-1,1],size(XxSims,1),1); %correctness in original context
    f2_sims = f2_sims(:);
    f3_sims = repmat([-ones(1,4),ones(1,4)],size(XxSims,1),1);
    f3_sims = f3_sims(:);
    rfx_sims = (1:size(XxSims,1))';
    f4_sims = repmat(rfx_sims,8,1);
    
    %%% regressions with  participant as random variable %%%
    formula = 'Y~Valence*Information*Correctness + (1 | Participant)';
    formulaPref = 'Y~Valence*Information*Correctness'; %no random intercept for preference as it has to be 50 anyway
    formula2 = 'Y~Valence*Information*Correctness + (1| Participant)';
    formulaSims = 'Y~Valence*Information*Correctness + (1 | Participant)';
    
    %%% tables for all 10 exps
    tbl_mod_all = table(XxMod(:),f1,f2,f3,f4,categorical(f5),'variableNames',vnames);
    tbl_data_all = table(Xx(:),f1,f2,f3,f4,categorical(f5),'variableNames',vnames);
     
    %%% use extended dataset for acc and pref
    if k_meas ==1
        lmAllModelExtended{k_meas} = fitlme(tbl_mod_all,formulaPref);
        lmAllExtended{k_meas} = fitlme(tbl_data_all,formulaPref);
    elseif k_meas <3
        lmAllModelExtended{k_meas} = fitlme(tbl_mod_all,formula);
        % data %
        lmAllExtended{k_meas} = fitlme(tbl_data_all,formula);
        lmAll2Extended{k_meas} = fitlme(tbl_data_all,formula2);
        FExtended{k_meas} = anova(lmAll2Extended{k_meas});
    end
    
    %%% tables using only conf exps%%%
    idcconfexp = find(f5<=5);
    % model %
    tbl_mod = table(XxMod(idcconfexp),f1(idcconfexp),f2(idcconfexp),f3(idcconfexp),f4(idcconfexp),categorical(f5(idcconfexp)),'variableNames',vnames);
    tbl_data= table(Xx(idcconfexp),f1(idcconfexp),f2(idcconfexp),f3(idcconfexp),f4(idcconfexp),categorical(f5(idcconfexp)),'variableNames',vnames);
    tbl_sims = table(XxSims(:),f1_sims,f2_sims,f3_sims,f4_sims,'variableNames',vnamesSims);
    
    if k_meas == 1
        lmAllModel{k_meas} = fitlme(tbl_mod,formulaPref);
        lmAll{k_meas} = fitlme(tbl_data,formulaPref);
        lmAllSims{k_meas} = fitlme(tbl_sims,formulaPref);

    else
        lmAllModel{k_meas} = fitlme(tbl_mod,formula);
        lmAll{k_meas} = fitlme(tbl_data,formula);
        lmAllSims{k_meas} = fitlme(tbl_sims,formulaSims);
    end
  
    %%% ANOVA intermediate options L25 G25
    
    if k_meas == 1
        Xx2 = Xx(:,[2,3,6,7]);
        f1 = repmat([1,2,1,2],size(Xx2,1),1); %val
        f1 = f1(:);
        f2 = repmat([1,1,2,2],size(Xx2,1),1); %information
        f2 = f2(:);
        rfx = (1:size(Xx2,1))';
        f4 = repmat(rfx,4,1);
        [~,tblInversion,stats]=anovan(Xx2(:),{f1,f2,f4},'model','interaction','random',3,'varnames',{'Valence','Information','Participant'},'Display', 'on');
        [cdata,m,h,nms] = multcompare(stats,'Dimension',[1,2],'Display','off');
        congroupletters = makeGroupsMultipleComparisons(cdata);
    end
    

    %%% individual regressions per experiment, with participants as random
    % variables %%%
    for iexp = [unique(f5)']
        if iexp>5 && k_meas>2
            continue
        end
        
        disp(['exp ',num2str(iexp)])
        tbl_mod_exp = tbl_mod_all(f5==iexp,:);
        tbl_data_exp = tbl_data_all(f5==iexp,:);
        if k_meas ==1
            formula = ['Y~Valence*Information*Correctness'];
        else
            formula = ['Y~Valence*Information*Correctness+ (1 |Participant)'];
        end
        lmAllExp{k_meas,iexp} = fitlme(tbl_data_exp,formula);
        lmAllExpModel{k_meas,iexp} = fitlme(tbl_mod_exp,formula);    
    end
end

%%% make matrices of effects per iv, dv, and exp
npar = numel(lmAll{1}.Coefficients.Estimate);
parsExp = nan(4,npar,12);
parsExpSE = nan(4,npar,12);
parsExpModel = nan(4,npar,12);
parsExpModelSE = nan(4,npar,12);
tExp = nan(4,npar,12);
tExpModel = nan(4,npar,12);

parsExpSims = nan(4,npar,12);
parsExpSimsSE = nan(4,npar,12);
tExpSims = nan(4,npar,12);

%%% change order of experiments: 1-5 no confidence, 6-10 cnfidence
expid = [6:10,1:5];

for imeas = 1:4
    for ipar = 1:npar
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
                parsExpSims(imeas,ipar,iexp) =stats.Estimate(ipar);
                parsExpSimsSE(imeas,ipar,iexp) = stats.SE(ipar);
                tExpSims(imeas,ipar,iexp) = stats.tStat(ipar);
                
            elseif (iexp ==12 && imeas<3)
                [b,~,stats] = fixedEffects(lmAllExtended{imeas},'DFMethod','Satterthwaite');
                parsExp(imeas,ipar,iexp) = stats.Estimate(ipar);
                parsExpSE(imeas,ipar,iexp) = stats.SE(ipar);
                tExp(imeas,ipar,iexp) = stats.tStat(ipar);
                
                [b,~,stats] = fixedEffects(lmAllModelExtended{imeas},'DFMethod','Satterthwaite');
                parsExpModel(imeas,ipar,iexp) = stats.Estimate(ipar);
                parsExpModelSE(imeas,ipar,iexp) = stats.SE(ipar);
                tExpModel(imeas,ipar,iexp) = stats.tStat(ipar);
                
                
            elseif ~(iexp>5 && imeas>2)
                thisexp = expid(iexp);
                [b,~,stats] = fixedEffects(lmAllExp{imeas,iexp},'DFMethod','Satterthwaite');
                parsExp(imeas,ipar,thisexp) = stats.Estimate(ipar);
                parsExpSE(imeas,ipar,thisexp) = stats.SE(ipar);
                tExp(imeas,ipar,thisexp) = stats.tStat(ipar);
                
                [b,~,stats] = fixedEffects(lmAllExpModel{imeas,iexp},'DFMethod','Satterthwaite');
                parsExpModel(imeas,ipar,thisexp) = stats.Estimate(ipar);
                parsExpModelSE(imeas,ipar,thisexp) = stats.SE(ipar);
                tExpModel(imeas,ipar,thisexp) = stats.tStat(ipar);
                
            end
        end
    end
end


%% Plot effects 

coldat = [.2,.2,.2];
colmod = [0.4,0.4,0.9];
colsims = [0.0,0.7,0.2];


%%% plots per measure, parameter and experiment

% main effects
figure()
measnames = {'Preference','Accuracy','Confidence','Overconfidence'};
parnames = {'B_{0}','B_{V}','B_{C}','B_{I}',...
    'B_{V:C}','B_{V:I}','B_{C:I}','B_{V:C:I}'};

for imeas = 1:4
    for ipar = 1:4
        subplot(4,4,sub2ind([4,4],ipar,imeas))
        hold on
        ylabel(parnames{ipar})
        title(['Transfer: ',measnames{imeas}]);
        l1 =errorbar([1:12]-0.2,squeeze(parsExp(imeas,ipar,:)),squeeze(parsExpSE(imeas,ipar,:)),'Marker','None','LineStyle','None','LineWidth',1.5,'Color',coldat);
        alpha(0.6)
        l2 = errorbar([1:12],squeeze(parsExpModel(imeas,ipar,:)),...
            squeeze(parsExpModelSE(imeas,ipar,:)),'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colmod);
        
        l3 = errorbar([1:12]+0.2,squeeze(parsExpSims(imeas,ipar,:)),...
            squeeze(parsExpSimsSE(imeas,ipar,:)),'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colsims);
        alpha(0.5)
        
        if numel(unique(sign(get(gca,'YLim'))))>1
            hl = plot([0,13],[0,0],'k--')
            uistack(hl,'bottom')
        end
        
        yl = get(gca,'YLim');
       
        mymax = max([squeeze(parsExp(imeas,ipar,:))+squeeze(parsExpSE(imeas,ipar,:));squeeze(parsExpModel(imeas,ipar,:))+squeeze(parsExpModelSE(imeas,ipar,:))]);
        mymin = min([squeeze(parsExp(imeas,ipar,:))+squeeze(parsExpSE(imeas,ipar,:));squeeze(parsExpModel(imeas,ipar,:))+squeeze(parsExpModelSE(imeas,ipar,:))]);
        
        for iexp = 1:12
            if iexp >=11
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
        
        vl = plot([10.5,10.5],yl,'k-');
        
        xticks(1:12)
        xticklabels({'1','2','3','4','5','6','7','8','9','10','Conf','All'});
        xtickangle(90)
        xlim([0,13])
        if imeas ==1 && ipar ==1
            leg = legend([l1,l2,l3],{'Data','Predicted','Simulation'},'Location','best');
            leg.BoxFace.ColorType='truecoloralpha';
            leg.BoxFace.ColorData=(uint8(255*[1 1 1 0.75]'));
        end
        set(gca,'FontName','Arial','FontSize',11);
    end
end
set(gcf,'Position',[300,100,1200,700])
saveas(gcf,['Plots/coeffMainTT',outsuffix,'.png']);
saveas(gcf,['Plots/coeffMainTT',outsuffix,'.svg']);

% interactions
figure()
for imeas = 1:4
    for ipar = 5:8
        subplot(4,4,sub2ind([4,4],ipar-4,imeas))
        hold on
        ylabel(parnames{ipar})
        l1 =errorbar([1:12]-0.2,squeeze(parsExp(imeas,ipar,:)),squeeze(parsExpSE(imeas,ipar,:)),'Marker','None','LineStyle','None','LineWidth',1.5,'Color',coldat);
        alpha(0.6)
        title(['Transfer :', measnames{imeas}]);
        l2 = errorbar([1:12],squeeze(parsExpModel(imeas,ipar,:)),...
            squeeze(parsExpModelSE(imeas,ipar,:)),'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colmod);
        
        l3 = errorbar([1:12]+0.2,squeeze(parsExpSims(imeas,ipar,:)),...
            squeeze(parsExpSimsSE(imeas,ipar,:)),'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colsims);
        if numel(unique(sign(get(gca,'YLim'))))>1
            hl = plot([0,13],[0,0],'k--')
            uistack(hl,'bottom')
        end
        
        mymax = max([squeeze(parsExp(imeas,ipar,:))+squeeze(parsExpSE(imeas,ipar,:));squeeze(parsExpModel(imeas,ipar,:))+squeeze(parsExpModelSE(imeas,ipar,:))]);
        mymin = min([squeeze(parsExp(imeas,ipar,:))+squeeze(parsExpSE(imeas,ipar,:));squeeze(parsExpModel(imeas,ipar,:))+squeeze(parsExpModelSE(imeas,ipar,:))]);
        
        yl = get(gca,'YLim');
        
        for iexp = 1:12
            if iexp >=11
                w = 'bold';
            else
                w = 'normal';
            end
            
            if abs(tExp(imeas,ipar,iexp))>1.96
                text(iexp,mymax+0.2*(yl(2)-yl(1)),'*','Color',coldat)
            end
            if abs(tExpModel(imeas,ipar,iexp))>1.96
                text(iexp,mymax+0.1*(yl(2)-yl(1)),'*','Color',colmod)
            end
        end
        
        if imeas ==1 && (ipar-4) ==1
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
        xtickangle(90)
        xlim([0,13])
        if imeas ==1 && (ipar-4) ==1
            leg = legend([l1,l2,l3],{'Data','Predicted','Simulation'},'Location','South');
            
            leg.BoxFace.ColorType='truecoloralpha';
            leg.BoxFace.ColorData=(uint8(255*[1 1 1 0.75]'));
        end
        set(gca,'FontName','Arial','FontSize',11)
    end
end
set(gcf,'Position',[300,100,1200,700])
saveas(gcf,['Plots/coeffInteractionTT',outsuffix,'.png']);
saveas(gcf,['Plots/coeffInteractionTT',outsuffix,'.svg']);


%% Plot sorted coefficients (pooling exps)

%%% all data, only acc/pref
figure()
% parnames = {'\beta_{0}','\beta_{Valence}','\beta_{Information}','\beta_{Interaction}'};
parnames = {'0','V','C','I',...
    'V:C','V:I','C:I','V:C:I'};

% xlabels = {'Data','\DeltaQ','\DeltaQ+Q_C','\DeltaQ+\SigmaQ','\DeltaQ+V'};
measnames = {'Preference','Accuracy','Confidence','Overconfidence'};
npar = numel(lmAll{1}.Coefficients.Estimate);

for imeas = 1:2
    subplot(4,1,imeas)
    hold on
    
    dumM = parsExp(imeas,2:end,12);
    [~,id] = sort((dumM),'descend');
    dumM = dumM(id);
    dumSE = parsExpSE(imeas,2:end,12);
    dumSE = dumSE(id);
    
    offset = 0.1;
    
    l1 =errorbar([1:npar-1]-offset,dumM,dumSE,'Marker','None','LineStyle','None','LineWidth',1.5,'Color',coldat);
    
    title(['Transfer: ',measnames{imeas},'(all exp.)']);
    ylabel('Regression fixed-effect');
    dummodelM = parsExpModel(imeas,2:end,12);
    dummodelSE = parsExpModelSE(imeas,2:end,12);
    dummodelM = squeeze(dummodelM(id));
    dummodelSE = squeeze(dummodelSE(id));
    
    l2 = errorbar([1:npar-1],dummodelM,dummodelSE,'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colmod);
    
    dumsimsM = parsExpSims(imeas,2:end,11);
    dumsimsSE = parsExpSimsSE(imeas,2:end,11);
    dumsimsM = squeeze(dumsimsM(id));
    dumsimsSE = squeeze(dumsimsSE(id));
    
    l3 = errorbar([1:npar-1]+offset,dumsimsM,dumsimsSE,'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colsims);
    
    if numel(unique(sign(get(gca,'YLim'))))>1
        hl = plot([0,9],[0,0],'k--')
        uistack(hl,'bottom')
    end
    
    yl = get(gca,'YLim');
    
    mymax = max(max([squeeze(parsExp(imeas,2:end,12))+squeeze(parsExpSE(imeas,2:end,12));squeeze(parsExpModel(imeas,2:end,12))'+squeeze(parsExpModelSE(imeas,2:end,12))]));
    mymin = min(min([squeeze(parsExp(imeas,2:end,12))+squeeze(parsExpSE(imeas,2:end,12));squeeze(parsExpModel(imeas,2:end,12))'-squeeze(parsExpModelSE(imeas,2:end,12))]));
    
    %show significance
    dumtExp = tExp(imeas,2:end,12);
    dumtExp = dumtExp(id);
    dumtExpModel = squeeze(tExpModel(imeas,2:end,12));
    dumtExpModel = dumtExpModel(id);
    
    for ipar = 1:numel(id)
        if abs(dumtExp(ipar))>1.96
            text(ipar,mymax+0.2*(yl(2)-yl(1)),'*','Color',coldat,'HorizontalAlignment','center')
        end
        
        if abs(dumtExpModel(ipar))>1.96
            text(ipar,mymax+0.1*(yl(2)-yl(1)),'*','Color',colmod,'HorizontalAlignment','center')
        end
    end
    
    yl = get(gca,'YLim');
    
    ylim([mymin,mymax+0.25*(diff(yl))])
    
    yl = get(gca,'YLim');
    if numel(unique(sign(get(gca,'YLim'))))>1
        plot([0,npar],[0,0],'k--')
        %         uistack(hl,'bottom')
    end
    
    yl = get(gca,'YLim');
    
    xticks(1:npar-1)
    dumxlabels = parnames(2:end);
    xlim([0,numel(dumxlabels)+1])
    
    xticklabels(dumxlabels(id))
    xtickangle(90)
    
    set(gca,'FontName','Arial','FontSize',11)
    
end

set(gcf,'Position',[0,0,400,1000])
saveas(gcf,'Plots/coefSortedTT_allData.svg');


%%% Confidence experiments only 
figure()

for imeas = 1:4
    subplot(4,1,imeas)
    hold on
    
    dumM = parsExp(imeas,2:end,11);
    [~,id] = sort((dumM),'descend');
    dumM = dumM(id);
    dumSE = parsExpSE(imeas,2:end,11);
    dumSE = dumSE(id);
    
    offset = 0.1;
    l1 =errorbar([1:npar-1]-offset,dumM,dumSE,'Marker','None','LineStyle','None','LineWidth',1.5,'Color',coldat);
    title(['Transfer: ',measnames{imeas}]);
    ylabel('Regression fixed-effect');
    dummodelM = parsExpModel(imeas,2:end,11);
    dummodelSE = parsExpModelSE(imeas,2:end,11);
    dummodelM = squeeze(dummodelM(id));
    dummodelSE = squeeze(dummodelSE(id));
    
    l2 = errorbar([1:npar-1],dummodelM,dummodelSE,'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colmod);
    
    dumsimsM = parsExpSims(imeas,2:end,11);
    dumsimsSE = parsExpSimsSE(imeas,2:end,11);
    dumsimsM = squeeze(dumsimsM(id));
    dumsimsSE = squeeze(dumsimsSE(id));
    
    l3 = errorbar([1:npar-1]+offset,dumsimsM,dumsimsSE,'Marker','None','LineStyle','None','LineWidth',1.5,'Color',colsims);
    
    if numel(unique(sign(get(gca,'YLim'))))>1
        hl = plot([0,9],[0,0],'k--')
        uistack(hl,'bottom')
    end
    
    yl = get(gca,'YLim');
    
    mymax = max(max([squeeze(parsExp(imeas,2:end,11))+squeeze(parsExpSE(imeas,2:end,11));squeeze(parsExpModel(imeas,2:end,11))'+squeeze(parsExpModelSE(imeas,2:end,11))]));
    mymin = min(min([squeeze(parsExp(imeas,2:end,11))+squeeze(parsExpSE(imeas,2:end,11));squeeze(parsExpModel(imeas,2:end,11))'-squeeze(parsExpModelSE(imeas,2:end,11))]));
    
    %show significance
    dumtExp = tExp(imeas,2:end,11);
    dumtExp = dumtExp(id);
    dumtExpModel = squeeze(tExpModel(imeas,2:end,11));
    dumtExpModel = dumtExpModel(id);
    
    for ipar = 1:numel(id)
        if abs(dumtExp(ipar))>1.96
            text(ipar-0,mymax+0.2*(yl(2)-yl(1)),'*','Color',coldat,'HorizontalAlignment','center')
        end
        
        if abs(dumtExpModel(ipar))>1.96
            text(ipar+0,mymax+0.1*(yl(2)-yl(1)),'*','Color',colmod,'HorizontalAlignment','center')
        end
    end
    
    yl = get(gca,'YLim');
    ylim([mymin,mymax+0.25*(diff(yl))])
    yl = get(gca,'YLim');
    if numel(unique(sign(get(gca,'YLim'))))>1
        plot([0,npar],[0,0],'k--')
    end
    
    yl = get(gca,'YLim');
    
    xticks(1:npar-1)
    dumxlabels = parnames(2:end);
    xlim([0,numel(dumxlabels)+1])
    
    xticklabels(dumxlabels(id))
    xtickangle(90) 
    set(gca,'FontName','Arial','FontSize',11)  
end

set(gcf,'Position',[0,0,400,1000])
saveas(gcf,'Plots/coefSortedTT_confData.svg');


%% Write latex output 
latexCode = [];

for imeas = 1:4
    
    %%%save tables in latex format
    %%% -iterate over measures
    %%% - save single tex file with all tables
    
    if imeas <3
        thisTable = convert_LMM2latex(...
            {lmAllExtended{imeas},lmAllModelExtended{imeas}},... % array of glme results
            {'Data','Model'},... % name of the models
            {'tStat'},... % information required
            2,... % round digits
            ['regTT',measnames{imeas},'_allexp'],... % table label
            'vertical'); % orientation (vertical or rotated)
        
        thisTable=replaceWords(thisTable,...
            {'Insert Title Here.'},... % old names
            {['Regression results for ',measnames{imeas}, ' in Transfer Task (all data)']}); % % new names
        
        latexCode =[latexCode,thisTable];
        
    end
    
    thisTable = convert_LMM2latex(...
        {lmAll{imeas},lmAllModel{imeas},lmAllSims{imeas}},... % array of glme results
        {'Data','Model','Simulations'},... % name of the models
        {'tStat'},... % information required
        2,... % round digits
        ['regTT',measnames{imeas},'_allexp'],... % table label
        'vertical'); % orientation (vertical or rotated)
    
    thisTable=replaceWords(thisTable,...
        {'Insert Title Here.'},... % old names
        {['Regression results for ',measnames{imeas}, ' in Transfer Task (confidence experiments)']}); % % new names
    
    latexCode =[latexCode,thisTable];
    
    thisTable = convert_LMM2latex(...
        lmAllExp(imeas,1:5),... % array of glme results
        {'Exp. 1', 'Exp. 2', 'Exp. 3', 'Exp. 4','Exp.5'},... % name of the models
        {'tStat'},... % information required
        2,... % round digits
        ['regTT',measnames{imeas},'_confexps'],... % table label
        'vertical'); % orientation (vertical or rotated)
    
    thisTable=replaceWords(thisTable,...
        {'Insert Title Here.'},... % old names
        {['Regression results for ',measnames{imeas}, ' in Transfer Task, for each confidence experiment']}); % % new names
    
    latexCode =[latexCode,thisTable];
    
    if imeas ==1
        thisTable = convert_LMM2latex(...
            lmAllExp(imeas,6:10),... % array of glme results
            {'Exp. 6', 'Exp. 7', 'Exp. 8', 'Exp. 9','Exp.10'},... % name of the models
            {'tStat'},... % information required
            2,... % round digits
            ['regTT',measnames{imeas},'_noconfexps'],... % table label
            'vertical'); % orientation (vertical or rotated)
        
        thisTable=replaceWords(thisTable,...
            {'Insert Title Here.'},... % old names
            {['Regression results for ',measnames{imeas}, ' in Transfer Task, for each non-confidence experiment']}); % % new names
        
        latexCode =[latexCode,thisTable];
        
    end
    
    thisTable = convert_LMM2latex(...
        lmAllExpModel(imeas,1:5),... % array of glme results
        {'Exp. 1', 'Exp. 2', 'Exp. 3', 'Exp. 4','Exp.5'},... % name of the models
        {'tStat'},... % information required
        2,... % round digits
        ['regTT',measnames{imeas},'Model_confexps'],... % table label
        'vertical'); % orientation (vertical or rotated)
    
    thisTable=replaceWords(thisTable,...
        {'Insert Title Here.'},... % old names
        {['Regression results for model-predicted ',measnames{imeas}, ' in Transfer Task, for each confidence experiment']}); % % new names
    
    latexCode =[latexCode,thisTable];
    
    if imeas ==1
        thisTable = convert_LMM2latex(...
            lmAllExpModel(imeas,6:10),... % array of glme results
            {'Exp. 6', 'Exp. 7', 'Exp. 8', 'Exp. 9','Exp.10'},... % name of the models
            {'tStat'},... % information required
            2,... % round digits
            ['regTT',measnames{imeas},'Model_noconfexps'],... % table label
            'vertical'); % orientation (vertical or rotated)
        
        thisTable=replaceWords(thisTable,...
            {'Insert Title Here.'},... % old names
            {['Regression results for model-predicted ',measnames{imeas}, ' in Transfer Task, for each non-confidence experiment']}); % % new names
        
        latexCode =[latexCode,thisTable];
        
    end
end

%%% write file
fid=fopen('Results/LatexTables/statsTT.tex','w'); fprintf(fid,'%s\n',latexCode{:}); fclose(fid);

