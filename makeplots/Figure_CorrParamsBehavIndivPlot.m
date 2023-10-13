%%% make scatterplot matrix of parameters (RL and confidence) against
%%% behavior (accuracy, confidence, valence bias, overconfidence),
%%% overlays lines from simulations of different models, and writes table
%%% of correlations

% clear

whichTasks = [1:2]; %1 Learning, 2 Transfer
whichPartype = [1,2]; %1 RL model, 2 Confidence model
whichPar = {[4,6],[1,3]}; %first vector for RL params, second confidence params
whichMeas = [1:4]; %1 acc, 2 conf, 3 valence bias, 4 overconf
whichModels = 2;

addpath('helperfuncs\')
addpath('helperfuncs\makeLatexTable\')

nlearnmodel = 11; %the learning model to use

%%% load parameters for all experiments
rl = load('Results\RLvars_all.mat');
% clear rl
params_learn = [cell2mat(rl.params_allexp(1:90,nlearnmodel))];

regLT = load('Results\reg_conflogit_learning_dqabs.mat','confcoeff','whichConfModel','whichLearnModel','isConfPrev','confBias')
regTT = load('Results\reg_conflogit_posttest_dqabs.mat','confcoeff','whichConfModel','whichLearnModel','isConfPrev','confBias')

idcLT = find(strcmp(regLT.confBias,'+Qc') & regLT.whichLearnModel == nlearnmodel & regLT.isConfPrev);
idcTT = find(strcmp(regTT.confBias,'+Qc')& regTT.whichLearnModel == nlearnmodel & ~regTT.isConfPrev);

confparamsLT =  regLT.confcoeff{idcLT};
confparamsTT =  regTT.confcoeff{idcTT};

load('Results/corrParamsBehav.mat'); %created by running corrParamsBehav.m
npreg = load('Results\NonParRegSims.mat');% load non-par reg
sims = load('Results\SimulationsSubjParams.mat');

%%% load data files
resultsdir = ['Results', filesep];
exps = load([resultsdir,'data_all'])
exps = exps.data_all;


%% load behavioral results
kk = 0;
% pref = [];
% accPost = [];
% confpost =[];
% conf_mat_sess = [];
% corr_mat_sess = [];

pref  = rl.pref(1:90,:);
[prefMat_all,accMat_all,confMat_all,confMatPair_all] = makePostMats(rl.ss,rl.aa,rl.confpost);
accPost = squeeze(nanmean(accMat_all(1:90,:,:)./4,3));
confpost = squeeze(nanmean(confMatPair_all(1:90,:,:),3));

conf_meanLT =  nanmean(nanmean(nanmean(rl.conf(1:90,:,:,:),2),3),4).*100;
acc_meanLT = nanmean(nanmean(nanmean(rl.correct(1:90,:,:,:),2),3),4).*100;
calib_meanLT = conf_meanLT-acc_meanLT;

yLT = squeeze(nanmean(nanmean(rl.conf(1:90,:,:,:).*100,2),4));
valbias_meanLT = nanmean(yLT(:,[1,2])-yLT(:,[3,4]),2);

conf_meanTT = nanmean(confpost,2).*100;
acc_meanTT = nanmean(accPost,2).*100;
calib_meanTT = conf_meanTT-acc_meanTT;
pair_mat = [6 5;2 1;8 7;4 3;];
yTT = [confpost(:,pair_mat(1,:)),confpost(:,pair_mat(2,:)),confpost(:,pair_mat(3,:)),confpost(:,pair_mat(4,:))];
valbias_meanTT = nanmean(yTT(:,[3,4,7,8])-yTT(:,[1,2,5,6]),2).*100;

allmeasLT = [acc_meanLT,conf_meanLT,valbias_meanLT,calib_meanLT];
allmeasTT = [acc_meanTT,conf_meanTT,valbias_meanTT,calib_meanTT];


%%% simulated measures 
% Learning
accSims = 100.*squeeze(nanmean(nanmean(sims.outSorted.pc,3),5));
confSims = 100.*squeeze(nanmean(nanmean(nanmean(sims.outConfLSorted(:,:,2,:,:,:),4),6),1));

acc_meanLT_sims = nanmean(nanmean(nanmean(accSims,2),3),4);
conf_meanLT_sims = nanmean(nanmean(nanmean(confSims,2),3),4);
valbias_meanLT_sims = nanmean(confSims(:,[1,2])-confSims(:,[3,4]),2);
calib_meanLT_sims = conf_meanLT_sims-acc_meanLT_sims;


% Transfer
accSimsT = 100.*squeeze(nanmean(sims.accPost));
confSimsT = 100.*squeeze(nanmean(sims.confPost(:,2,:,:)));
   
conf_meanTT_sims = nanmean(confSimsT,2);
acc_meanTT_sims = nanmean(accSimsT,2);
calib_meanTT_sims = conf_meanTT_sims-acc_meanTT_sims;
pair_mat = [6 5;2 1;8 7;4 3];
newIdc = pair_mat';
yTT = [confSimsT(:,pair_mat(1,:)),confSimsT(:,pair_mat(2,:)),confSimsT(:,pair_mat(3,:)),confSimsT(:,pair_mat(4,:))];
valbias_meanTT_sims = nanmean(yTT(:,[3,4,7,8])-yTT(:,[1,2,5,6]),2).*100;

allmeasLT_sims = [acc_meanLT_sims,conf_meanLT_sims,valbias_meanLT_sims,calib_meanLT_sims];
allmeasTT_sims = [acc_meanTT_sims,conf_meanTT_sims,valbias_meanTT_sims,calib_meanTT_sims];


% accPost_sim = squeeze(nanmean(sims.accPost,2)

bias = (params_learn(:,2)-params_learn(:,3))./(params_learn(:,2)+params_learn(:,3));
bias2 = (params_learn(:,2)-params_learn(:,3));%non-normalized

learnparams_bias = [params_learn,bias,bias2];

bias_sims =  (sims.learnparams(:,2)-sims.learnparams(:,3))./(sims.learnparams(:,2)+sims.learnparams(:,3));
bias2_sims =  (sims.learnparams(:,2)-sims.learnparams(:,3));
learnparams_bias_sims = [sims.learnparams,bias_sims,bias2_sims];

modelcolor = [102,194,165;
    252,141,98;
    141,160,203]/255;

tasknames = {'Learning','Transfer'};
paramtypesnames = {'RL','conf'};
measnames = {'Accuracy','Confidence','Valence bias','Overconfidence'};
learnparamnamesTable = {'$log(\beta)$','$\alpha_{CON}$','$\alpha_{DIS}$','$log(\alpha_V)$','w','$\frac{\alpha_{CON}-\alpha_{DIS}}{\alpha_{CON}+\alpha_{DIS}}$','$\alpha_{CON}-\alpha_{DIS}$'};
learnparamnames = {'logTemperature','alpaCON','alphaDIS','alphaV','w','lbias','lbiascontrast'};
confparamnamesTable = {'\beta_0','$\beta_{|\Delta Q|}$','$\beta_{BIAS}$','$\beta_{y(t-1)}$'};
confparamnames = {'b0','bDQ','bBias','bConfPrev'};
learnparamnamesSaving = {'beta','alpaCON','alphaDIS','alphaV','w','lbias','lbiascontrast'};
confparamnamesSaving = {'b0','bDQ','bBias','bConfPrev'};
confmodellabels = {'Q_C','\SigmaQ','V'};

for itask = whichTasks;
    switch itask
        case 1
            allmeas = allmeasLT;
            myT = 'Learning';
            confparams = confparamsLT;
            nmodels = numel(idcLT);
        case 2
            allmeas = allmeasTT;
            myT = 'Transfer';
            confparams = confparamsTT;
            nmodels = numel(idcTT);
    end
    for ipartype = whichPartype
        switch ipartype
            case 1
                params = learnparams_bias;
                params(:,1) = log(params(:,1));
                paramnames = [learnparamnames];
                paramnamesSaving = learnparamnamesSaving;
            case 2
                params = confparams;
                paramnames = confparamnames;
                paramnamesSaving = confparamnamesSaving;
        end
        
        for ipar = whichPar{ipartype}
            X = params(:,ipar);
            Xreg = X;
        	Xline = linspace(min(X),max(X),1000)';
            XlineReg = Xline;

            if ipartype ==1 && ipar ==4
                Xreg = log(X);
%                 Xline = log(Xline);
                XlineReg = log(XlineReg);

            end
            
            for imeas = whichMeas
                Y = allmeas(:,imeas);
                figure()
                set(gcf,'Color',[1,1,1])
                subplot('position',[0.82,0.15,0.12,0.6])
                h = histogram(Y,15,'Orientation','horizontal','FaceColor',0.7*[1,1,1],'EdgeColor',[1,1,1],...
                    'BinLimits',[min(Y),max(Y)]);
                ylim([min(Y),max(Y)])
                yticks({})
                box off
                subplot('Position',[0.15,0.15,0.6,0.6])
                hold on
                % ylim([min([X;Y]),max([X;Y])])
                % xlim([min([X;Y]),max([X;Y])])
                xlim([min(X),max(X)])
                ylim([min(Y),max(Y)])
                % plot([min(X),max(X)],[min(Y),max(Y)],'Color',[.5,.5,.5],'LineStyle',':')
                
                %%% plot simulation nonpar regression line
                colnpreg = [102,194,165]/255;
                coldatafit = [252,141,98]/255;
                colsimfit =  [141,160,203]/255;
                residCV = cell(3,2);
                residLinCV = cell(3,2);
                legendR2CV = confmodellabels;
                for iconfmodel = 1:numel(whichModels)
                    thisModel = whichModels(iconfmodel);
                    dumX = npreg.regAllX{thisModel,itask,ipartype}{ipar}';
                    f_all = npreg.regAllF{thisModel,itask,ipartype}{ipar,imeas};
                    f = nan(size(Y));
                    for iY = 1:numel(Y)
                        [Xdist{itask,ipartype}{ipar,imeas}(iY),idc] = min(abs(dumX-X(iY)));
                        f(iY) = f_all(idc);
                    end
                    
                    sel = find(X<=max(dumX) & X>=min(dumX))%select points within range of simulation npreg
%                     sel = 1:numel(X);
                    residCV{iconfmodel,itask} = Y(sel)-f(sel);
                    
                    ssRes = sum(residCV{iconfmodel,itask}.^2);
                    ssReg = sum((f(sel)-mean(Y(sel))).^2);
                    ssTot =  sum((Y(sel)-mean(Y(sel))).^2);
                    
                    
                    fLin = npreg.regLin{thisModel,itask,ipartype}{ipar,imeas}.predict(Xreg);
                    residLinCV{iconfmodel,itask} = Y-fLin;
                    ssResLin = sum(residLinCV{iconfmodel,itask}.^2);
                    ssRegLin = sum((fLin-mean(Y)).^2);
                    ssTotLin =  sum((Y-mean(Y)).^2);
                    
                    
                    yci = npreg.regAllCI95{thisModel,itask,ipartype}{ipar,imeas};
                    l(iconfmodel) = plot(dumX,f_all,'Color',colnpreg);
                    fill([dumX;flipud(dumX)],[yci(2,:)';flipud(yci(1,:)')],colnpreg,'LineStyle','None','FaceAlpha',0.2);
                    rmse(iconfmodel,itask) = sqrt(mean(residCV{iconfmodel,itask}.^2));
                    rsquared(iconfmodel,itask) =  1- ssRes./ssTot;
                    
                    rsquaredLin(iconfmodel,itask) =  1- ssResLin./ssTotLin;
                    
                    %                     varExp=  (ssReg)/((ssRes+ssReg));
                    legendData = sprintf('data (\\rho = %0.2f, p = %0.2f)',rDataAll{itask,ipartype}(ipar,imeas),pDataAll{itask,ipartype}(ipar,imeas));
                    legendSim = sprintf('sim (\\rho = %0.2f, p = %0.2f)',npreg.rho{thisModel,itask,ipartype}(ipar,imeas),npreg.p{thisModel,itask,ipartype}(ipar,imeas));
                    legendR2CV = sprintf('CV non-param: R^2 = %0.2f',rsquared(iconfmodel,itask)),')';
%                     legendR2LinCV = sprintf('CV param.: R^2 = %0.2f',rsquaredLin(iconfmodel,itask)),')';

                end
                datapoints = scatter(X,Y,'Marker','o','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1]);
                
%                 dumX = [min([X;Y]):0.1:max([X;Y])]';
                
                thislm = fitlm(Xreg,Y,'RobustOpts','on');
                lmAllData{thisModel,itask,ipartype}{ipar,imeas} = thislm;
                [y,yci] = thislm.predict(XlineReg);    
                
%                 if ipar ==4 && ipartype==1
%                     [ySim,yciSim] = npreg.regLin{thisModel,itask,ipartype}{ipar,imeas}.predict(Xline);    
%                 else
                [ySim,yciSim] = npreg.regLin{thisModel,itask,ipartype}{ipar,imeas}.predict(XlineReg);    
%                 end
                
                
                r2data = thislm.Rsquared.Ordinary;
                legendR2data = sprintf('data fit (R^2 = %0.2f)',r2data);
                legendR2CVLin = sprintf('sim fit lm (CV: R^2 = %0.2f)',rsquaredLin(iconfmodel,itask));

                
                ldata = plot(Xline,y,'Color',coldatafit);
                
                lSim = plot(Xline,ySim,'Color',colsimfit);

                 %     plot(dumX,yci,'LineStyle','--','Color',[0.5,0.0,0.5]);
                patch([Xline;flipud(Xline)],[yci(:,2);flipud(yci(:,1))],coldatafit,'LineStyle','None','FaceAlpha',0.2);

                patch([Xline;flipud(Xline)],[yciSim(:,2);flipud(yciSim(:,1))],colsimfit,'LineStyle','None','FaceAlpha',0.2);

                %                 leg = legend(l,thisLegend{1:nmodels},'location','best');
                ldum =  plot(0,0,'Color','white')
                leg = legend([datapoints,ldata,l,ldum,lSim],...
                    {legendData,legendR2data,legendSim,legendR2CV,legendR2CVLin},'location','best');
                
%                 legend('boxoff')
%                 title(leg,'Simulations')
                ylabel(measnames{imeas},'FontSize',16)
                xlabel(paramnames{ipar},'FontSize',16,'Interpreter','latex')
                set(gca,'FontSize',16)
                subplot('Position',[0.15,0.82,0.6,0.12])
                histogram(X,15,'FaceColor',0.7*[1,1,1],'EdgeColor',[1,1,1],...
                    'BinLimits',[min(X),max(X)]);
                xlim([min(X),max(X)])
                xticks({})
                box off
                set(gcf,'Position',[100,100,650,650]);
                saveas(gcf,['Plots\regNonPar',...
                    paramnamesSaving{ipar},measnames{imeas},tasknames{itask},'.svg'])
% %                 close()
            end
        end
    end
end


%%% multiple regression %%%

params_learn_meas = params_learn;
params_learn_meas(:,4) = log(params_learn_meas(:,4));
params_learn_meas(:,1) = log(params_learn_meas(:,1));
params_learn_meas_sims = sims.learnparams;
params_learn_meas_sims(:,4) = log(params_learn_meas_sims(:,4));
params_learn_meas_sims(:,1) = log(params_learn_meas_sims(:,1));

for itask = 1:2
    switch itask
        case 1
            allmeas = allmeasLT;
            confparams = confparamsLT;
            allmeas_sims = allmeasLT_sims;
            confparams_sims = sims.confparamsLT{2};
        case 2
            allmeas = allmeasTT;
            confparams = confparamsTT;
            allmeas_sims = allmeasTT_sims;
            confparams_sims = sims.confparamsTT{2};
    end
    
     for ipartype = whichPartype
        
        for imeas = 1:4
            Y = allmeas(:,imeas);
            Ysim = allmeas_sims(:,imeas);
            
            switch ipartype
                case 1
                    params = params_learn_meas;
                    params_sim = params_learn_meas_sims;
                    paramnames = learnparamnames(1:5);
                case 2
                    params = confparams;%(:,1:3);
                    params_sim = confparams_sims;%(:,1:3);
                    paramnames = confparamnames(1:4);
            end

            
                   
            varnames = {'Y',paramnames{1:size(params,2)}};
            formula = [varnames{1},'~ 1'];
            for ivar = 2:numel(varnames)
                formula = [formula,'+',varnames{ivar}];
            end

           tbl = array2table([Y,params],...
            'VariableNames',varnames);
            thislm = fitlm(tbl,formula,'RobustOpts','on');
%             thislm = fitlm(params,Y,'RobustOpts','on');
            lmAll{itask,ipartype,imeas} = thislm;
            coeffs{itask,ipartype,imeas} = thislm.Coefficients.Estimate;
            tbl = array2table([Y./std(Y),params./std(params)],...
            'VariableNames',varnames);
            thislm_stand = fitlm(tbl,formula,'RobustOpts','on');
            lmAllStand{itask,ipartype,imeas} = thislm_stand;
            coeffs_se{itask,ipartype,imeas} = thislm.Coefficients.SE;
            coeffs_stand{itask,ipartype,imeas} =thislm_stand.Coefficients.Estimate;
            coeffs_stand_se{itask,ipartype,imeas} =thislm_stand.Coefficients.SE;
            
            %%% regression on sims
            tbl = array2table([Ysim,params_sim],...
            'VariableNames',varnames);
            thislm = fitlm(tbl,formula,'RobustOpts','on');
%             thislm = fitlm(params_sim,Ysim,'RobustOpts','on');
            lmAll_sims{itask,ipartype,imeas} = thislm;
            coefs_sims{itask,ipartype,imeas} = thislm.Coefficients.Estimate;
        
            tbl = array2table([Ysim./std(Ysim),params_sim./std(params_sim)],...
            'VariableNames',varnames);
            thislm_stand_sims = fitlm(tbl,formula,'RobustOpts','on');
            lmAllStand_sims{itask,ipartype,imeas} = thislm_stand_sims;
            coeffs_se_sims{itask,ipartype,imeas} = thislm.Coefficients.SE;
            coeffs_stand_sims{itask,ipartype,imeas} =thislm_stand_sims.Coefficients.Estimate;
            coeffs_stand_se_sims{itask,ipartype,imeas} =thislm_stand_sims.Coefficients.SE;
        end
     end
end


%%make latex table for multiple regression
%%% columns for LT, TT, LT sim, TT sim
%%% one table per measure (val bias and overconf) and partype (learn, conf)

latexCode = '';
fid=fopen('Results/LatexTables/stats_multipleRegBiases.tex','w');
partypenames = {'Learning','Confidence'};
coeftypenames = {'Raw','Standardized'};

for icoeftype = 2 %raw vs standardized
    for ipartype = 1:2
        for imeas = 1:4
            if icoeftype==1
            dum ={lmAll{1:2,ipartype,imeas},lmAll_sims{1:2,ipartype,imeas}};
            else
            dum ={lmAllStand{1:2,ipartype,imeas},lmAllStand_sims{1:2,ipartype,imeas}};
            end
            thisTable = convert_regression2latex(dum,...
                {'Learning', 'Transfer', 'Learning sim','Transfer sim'},... % name of the models
                {'tStat'},... % information required
                2,... % round digits
                ['Multiple regression',measnames{imeas}],... % table label
                'vertical'); % orientation (vertical or rotated)
            %         latexCode =strcat(latexCode,thisTable{:});
            
                thisTable=replaceWords(thisTable,...
                {'Insert Title Here.'},... % old names
                {['Multiple regression of ',measnames{imeas},' on ',partypenames{ipartype},' parameters (', coeftypenames{icoeftype},' coefficients)']}); % % new names
                thisTable = replaceWords(thisTable,learnparamnames,learnparamnamesTable);
                thisTable = replaceWords(thisTable,confparamnames,confparamnamesTable);

            fprintf(fid,'%s\n',thisTable{:});
            
        end
    end
end
fclose(fid);

% fid=fopen('Results/LatexTables/stats_multipleRegBiases.tex','w'); fprintf(fid,'%s\n',latexCode{:}); fclose(fid);

learnparamnames = {'logTemperature','alpaCON','alphaDIS','alphaV','w','lbias','lbiascontrast'};
partypenames = {'learn','conf'};
for ipartype = 1:2
    figure()
    switch ipartype
        case 1
            paramnames = [learnparamnames];
            paramnamesSaving = learnparamnamesSaving;
        case 2
            paramnames = [confparamnamesSaving];
            paramnamesSaving = confparamnamesSaving;
    end
    
    for imeas = 3:4
        %     theseCoeffs = cell2mat(coeffs(:,imeas)');
        %     theseCoeffs_se = cell2mat(coeffs_se(:,imeas)');
        
%         if imeas ==3 & ipartype==1
%             paramnames = {'log(\beta)','\alpha_{CON}','\alpha_{DIS}','log(\alpha_V)','w'};
%         elseif imeas ==4 & ipartype ==1
%             paramnames = {'\beta','\alpha_{CON}','\alpha_{DIS}','\alpha_V','w'};
%         end
        
        subplot(1,2,imeas-2)
        title(measnames{imeas})
        if imeas == 3
            ylabel('Standardized coefficient')
        end
        xtickangle(90)
        xlim([0,numel(paramnames)+1])
        %         ylim([-100,100])
        hold on

        ncoefL = numel(coeffs_stand{1,ipartype,imeas}(2:end));
        ncoefT = numel(coeffs_stand{2,ipartype,imeas}(2:end));
        
        bars1 = bar([1:ncoefL]-0.15,coeffs_stand{1,ipartype,imeas}(2:end),'FaceColor',.4*[1,1,1],'BarWidth',0.3)
        bars2 = bar([1:ncoefT]+0.15,coeffs_stand{2,ipartype,imeas}(2:end),'FaceColor',.7*[1,1,1],'BarWidth',0.3)
        
        errorbar([1:ncoefL]-0.15,coeffs_stand{1,ipartype,imeas}(2:end),coeffs_stand_se{1,ipartype,imeas}(2:end),'k','LineStyle','none','CapSize',1)
        errorbar([1:ncoefT]+0.15,coeffs_stand{2,ipartype,imeas}(2:end),coeffs_stand_se{2,ipartype,imeas}(2:end),'k','LineStyle','none','CapSize',1)
        
        bars3 = errorbar([1:ncoefL]-0.15,coeffs_stand_sims{1,ipartype,imeas}(2:end),coeffs_stand_se_sims{1,ipartype,imeas}(2:end),...
        'Color',[0.5,0,0.7],'MarkerEdgeColor',[0.5,0,0.7],'MarkerFaceColor',.4*[1,1,1],'Marker','o','MarkerSize',4,'LineStyle','none','CapSize',5)
        errorbar([1:ncoefT]+0.15,coeffs_stand_sims{2,ipartype,imeas}(2:end),coeffs_stand_se_sims{2,ipartype,imeas}(2:end),...
        'Color',[0.5,0,0.7],'MarkerEdgeColor',[0.5,0,0.7],'MarkerFaceColor',.7*[1,1,1],'Marker','o','MarkerSize',4,'LineStyle','none','CapSize',5)                
%         scatter([1:ncoefL]-0.15,coeffs_stand_sims{1,ipartype,imeas}(2:end),'ro')
%         scatter([1:ncoefT]+0.15,coeffs_stand_sims{2,ipartype,imeas}(2:end),'ro')      

        set(gca,'Xtick',1:ncoefL,...
            'XTickLabel',paramnames(1:ncoefL),'FontSize',10,'FontName','Arial')
        %                 set(gca,'TickLabelInterpreter','latex')
        
        if imeas ==4
            legend([bars1,bars2,bars3],{'Learning','Transfer','Simulation'},'Location','best');
        end
    end
    
    set(gcf,'Position',[0,0,800,300])
    saveas(gcf,['Plots/corrBehavParamsBars_',partypenames{ipartype},'.svg']);
end


for ipartype = 1:2
    figure()
    switch ipartype
        case 1
            paramnames = [learnparamnames];
            paramnamesSaving = learnparamnamesSaving;
        case 2
            paramnames = [confparamnamesSaving];
            paramnamesSaving = confparamnamesSaving;
    end
    
    for imeas = 1:2
        
        subplot(1,2,imeas)
        title(measnames{imeas})
        if imeas == 1
            ylabel('Standardized coefficient')
        end
        xtickangle(90)
        xlim([0,numel(paramnames)+1])
        %         ylim([-100,100])
        hold on

        ncoefL = numel(coeffs_stand{1,ipartype,imeas}(2:end));
        ncoefT = numel(coeffs_stand{2,ipartype,imeas}(2:end));
        
        bars1 = bar([1:ncoefL]-0.15,coeffs_stand{1,ipartype,imeas}(2:end),'FaceColor',.4*[1,1,1],'BarWidth',0.3)
        bars2 = bar([1:ncoefT]+0.15,coeffs_stand{2,ipartype,imeas}(2:end),'FaceColor',.7*[1,1,1],'BarWidth',0.3)
        
        errorbar([1:ncoefL]-0.15,coeffs_stand{1,ipartype,imeas}(2:end),coeffs_stand_se{1,ipartype,imeas}(2:end),'k','LineStyle','none','CapSize',1)
        errorbar([1:ncoefT]+0.15,coeffs_stand{2,ipartype,imeas}(2:end),coeffs_stand_se{2,ipartype,imeas}(2:end),'k','LineStyle','none','CapSize',1)
        
        bars3 = errorbar([1:ncoefL]-0.15,coeffs_stand_sims{1,ipartype,imeas}(2:end),coeffs_stand_se_sims{1,ipartype,imeas}(2:end),...
        'Color',[0.5,0,0.7],'MarkerEdgeColor',[0.5,0,0.7],'MarkerFaceColor',.4*[1,1,1],'Marker','o','MarkerSize',4,'LineStyle','none','CapSize',5)
        errorbar([1:ncoefT]+0.15,coeffs_stand_sims{2,ipartype,imeas}(2:end),coeffs_stand_se_sims{2,ipartype,imeas}(2:end),...
        'Color',[0.5,0,0.7],'MarkerEdgeColor',[0.5,0,0.7],'MarkerFaceColor',.7*[1,1,1],'Marker','o','MarkerSize',4,'LineStyle','none','CapSize',5)                
%         scatter([1:ncoefL]-0.15,coeffs_stand_sims{1,ipartype,imeas}(2:end),'ro')
%         scatter([1:ncoefT]+0.15,coeffs_stand_sims{2,ipartype,imeas}(2:end),'ro')      

        set(gca,'Xtick',1:ncoefL,...
            'XTickLabel',paramnames(1:ncoefL),'FontSize',10,'FontName','Arial')
        %                 set(gca,'TickLabelInterpreter','latex')
        
        if imeas ==4
            legend([bars1,bars2,bars3],{'Learning','Transfer','Simulation'},'Location','best');
        end
    end
    
    set(gcf,'Position',[0,0,800,300])
    saveas(gcf,['Plots/corrBehavParamsBars2_',partypenames{ipartype},'.svg']);
end
