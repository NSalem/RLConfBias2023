% plot correlation matrix of learning rates between
% models.
% Matrices 1-3: 
% confirmatory, disconfirmatory, confirmatory - disconfirmatory. 
% diagonal shows correlation with overconfidence
% Matrix 4: 
% contextual learning rate, diagonal shows correlation with valence bias.

addpath('helperfuncs');
load('Results\RLvars.mat')
load('Results\data_all_CONF.mat');

%% load data
whichexp = 1:5;
conf_mat_sess = [];
prevconf_mat_sess = [];
corr_mat_sess = [];
pref = [];
accPost = [];
confpost =[];
prefMat_all= [];
confMat_all = [];
accMat_all = [];
for iexp =whichexp
    whichsess = 1:3;
    subjects    = data_all(iexp).subjects;           % 3 sess
    n_sub       = length(subjects);
    n_sess      = numel(whichsess);%size(data_all(iexp).con,2);
    n_trials    = size([data_all(iexp).con],3);
    expN        = data_all(iexp).expN;
    
    thisExp =  MakeMatrix_sess_cond_trl(data_all(iexp),1);
    conf_mat_sess = [conf_mat_sess;thisExp.conf_mat_sess.*100];
    prevconf_mat_sess = [prevconf_mat_sess;thisExp.prevconf_mat_sess.*100];
    corr_mat_sess = [corr_mat_sess;thisExp.corr_mat_sess];
    pref = [pref;thisExp.pref];
    confpost = [confpost;thisExp.confpostcon];
    accPost = [accPost;thisExp.accPost];
end

overconf = squeeze(nanmean(nanmean(nanmean(conf_mat_sess-corr_mat_sess*100,4),3),2));

%% corr alphaCON, alphaDIS, and difference. Diagonal shows corr w/overconf

idc = [8:12];
params_mymodels = params_allexp(:,idc);
myCorrMatLR1 = nan([numel(idc),numel(idc)]);
myCorrMatLR2 = nan([numel(idc),numel(idc)]);
myCorrMatLRDiff = nan([numel(idc),numel(idc)]);
pCorrMatLR1 = nan([numel(idc),numel(idc)]);
pCorrMatLR2 = nan([numel(idc),numel(idc)]);
pCorrMatLRDiff = nan([numel(idc),numel(idc)]);

lrModel1 = nan([size(params_mymodels),2]);
lrModel2 = lrModel1;
for imodel1 = 1:size(params_mymodels,2)
    for imodel2 = 1:size(params_mymodels,2)
        clear lrModel1 lrModel2
        if imodel1 == imodel2
            for isub = 1:size(params_mymodels,1)
                lrModel1(isub,:) = params_mymodels{isub,imodel1}([2:3]);
            end
            [myCorrMatLR1(imodel1,imodel2),pCorrMatLR1(imodel1,imodel2)] = corr(lrModel1(:,1),overconf);
            [myCorrMatLR2(imodel1,imodel2),pCorrMatLR2(imodel1,imodel2)] = corr(lrModel1(:,2),overconf);
            [myCorrMatLRDiff(imodel1,imodel2),pCorrMatLRDiff(imodel1,imodel2)] = corr(lrModel1(:,1)-lrModel1(:,2),overconf);

        else
            for isub = 1:size(params_mymodels,1)
                lrModel1(isub,:) = params_mymodels{isub,imodel1}([2:3]);
                lrModel2(isub,:) = params_mymodels{isub,imodel2}([2:3]);
            end
            [myCorrMatLR1(imodel1,imodel2),pCorrMatLR1(imodel1,imodel2)] = corr(lrModel1(:,1),lrModel2(:,1));
            [myCorrMatLR2(imodel1,imodel2),pCorrMatLR2(imodel1,imodel2)] = corr(lrModel1(:,2),lrModel2(:,2));
            [myCorrMatLRDiff(imodel1,imodel2),pCorrMatLRDiff(imodel1,imodel2)] = corr(lrModel1(:,1)-lrModel1(:,2),lrModel2(:,1)-lrModel2(:,2));
        end
    end
end

%%% plots
figure()
modelnames = {'CONF','RELCONF','RELCONF_{Q_U}','RELCONF_{R_{OTHER}}','RELCONF_{R_{LAST}}'}

for iplot = 1:3
    subplot(1,3,iplot)
    if iplot ==1
        X = myCorrMatLR1;
        P = pCorrMatLR1;
        myTitle = '\alpha_{CON}';
    elseif iplot == 2
        X = myCorrMatLR2;
        P = pCorrMatLR2;
        myTitle= '\alpha_{DIS}';
    elseif iplot ==3
        myTitle='\alpha_{CON}-\alpha_{DIS}';
        X = myCorrMatLRDiff;
        P = pCorrMatLRDiff;
    end
    X = tril(X);
    X(X==0) = NaN;
    b = imagesc(X);
    set(b,'AlphaData',~isnan(X));
    title(myTitle);
    for imodel1 = 1:size(X,1)
        for imodel2 = 1:size(X,2)
            
            if P(imodel1,imodel2)<0.001
                sigSymb = '***';
            elseif P(imodel1,imodel2)<=0.01
                sigSymb = '**';
            elseif P(imodel1,imodel2) <0.05
                sigSymb = '*';
            elseif P(imodel1,imodel2) >0.05
                sigSymb = '';
            end
            
            if imodel2>imodel1
                text(imodel1,imodel2,[num2str(round(X(imodel2,imodel1),2)),sigSymb],'HorizontalAlignment','center');
            elseif imodel1 == imodel2
                text(imodel1,imodel2,[num2str(round(X(imodel2,imodel1),2)),sigSymb],'Color','k','HorizontalAlignment','center');
            end
        end
    end
    
    colormap('jet');
    colorbar();
    caxis([-1,1]);
    %%% add legend
    xticklabels(modelnames);
    xtickangle(90);
    yticklabels(modelnames);
    xticks(1:numel(modelnames))
    yticks(1:numel(modelnames))
    %%% add text
end


set(gcf,'Position',[0,0,1800,400])
% saveas(gcf,'Plots\corrLearningRatesAcrossModels.png');
% saveas(gcf,'Plots\corrLearningRatesAcrossModels.svg');

%% corr alphaV, diagonal shows corr with valence bias in confidence
confbias = squeeze(nanmean(nanmean(nanmean(conf_mat_sess(:,:,[1,2],:)-conf_mat_sess(:,:,[3,4],:),2),3),4));
myCorrMatLR3 = nan([numel(idc),numel(idc)]);
pCorrMatLR3 = nan([numel(idc),numel(idc)]);

idc = [2:5,9:12];
params_mymodels = params_allexp(:,idc);
lrModel1 = nan([size(params_mymodels),2]);
lrModel2 = lrModel1;
for imodel1 = 1:size(params_mymodels,2)
    for imodel2 = 1:size(params_mymodels,2)
        clear lrModel1 lrModel2
        if imodel1 == imodel2
            for isub = 1:size(params_mymodels,1)
                lrModel1(isub) = params_mymodels{isub,imodel1}([4]);
            end
            [myCorrMatLR3(imodel1,imodel2),pCorrMatLR3(imodel1,imodel2)] = corr(lrModel1',confbias);
        else
            for isub = 1:size(params_mymodels,1)
                lrModel1(isub,:) = params_mymodels{isub,imodel1}([4]);
                lrModel2(isub,:) = params_mymodels{isub,imodel2}([4]);
            end
            [myCorrMatLR3(imodel1,imodel2),pCorrMatLR3(imodel1,imodel2)] = corr(lrModel1,lrModel2);
        end
    end
end



modelnames = {'REL','REL_{Q_U}','REL_{OTHER}','REL_{LAST}','CONFREL','CONFREL_{Q_U}','CONFREL_{OTHER}','CONFREL_{LAST}'}
figure()
X = myCorrMatLR3;
P = pCorrMatLR3;

X = tril(X)
X(X==0) = NaN;
b = imagesc(X);
set(b,'AlphaData',~isnan(X));
title('\alpha_V');
for imodel1 = 1:size(X,1)
    for imodel2 = 1:size(X,2)
        
        if P(imodel1,imodel2)<0.001
            sigSymb = '***';
        elseif P(imodel1,imodel2)<=0.01
            sigSymb = '**';
        elseif P(imodel1,imodel2) <0.05
            sigSymb = '*';
        elseif P(imodel1,imodel2) >0.05
            sigSymb = '';
        end
        
        if imodel2>imodel1
            text(imodel1,imodel2,[num2str(round(X(imodel2,imodel1),2)),sigSymb],'HorizontalAlignment','center');
        elseif imodel1 == imodel2
            text(imodel1,imodel2,[num2str(round(X(imodel2,imodel1),2)),sigSymb],'Color','k','HorizontalAlignment','center');
        end
    end
end

colormap('jet');
colorbar();
caxis([-1,1]);
%%% add legend
xticklabels(modelnames);
xtickangle(90);
yticklabels(modelnames);
xticks(1:numel(modelnames))
yticks(1:numel(modelnames))

% saveas(gcf,'Plots\corrAlphaVAcrossModels.png');
% saveas(gcf,'Plots\corrAlphaVAcrossModels.svg');

