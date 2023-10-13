%%% make scatterplot matrix of parameters (RL and confidence) against
%%% behavior (accuracy, confidence, valence bias, overconfidence),
%%% overlays lines from simulations of different models, and writes table
%%% of correlations

clear

addpath('helperfuncs\')

nlearnmodel = 11; %the learning model to use

rlVars = load('Results\RLvars_all.mat');
corr_mat_sess = rlVars.correct;
conf_mat_sess = rlVars.conf;


regLT = load('Results\reg_conflogit_learning_dqabs.mat','confcoeff','whichConfModel','whichLearnModel','isConfPrev','confBias');
regTT = load('Results\reg_conflogit_posttest_dqabs.mat','confcoeff','whichConfModel','whichLearnModel','isConfPrev','confBias');
% idcLT = find(regLT.whichConfModel == 2 & regLT.whichLearnModel == nlearnmodel & regLT.isConfPrev);
% idcTT = find(regTT.whichConfModel == 2 & regTT.whichLearnModel == nlearnmodel & ~regTT.isConfPrev);

idcLT = find(strcmp(regLT.confBias,'+Qc') & regLT.whichLearnModel == nlearnmodel & regLT.isConfPrev);
idcTT = find(strcmp(regTT.confBias,'+Qc')& regTT.whichLearnModel == nlearnmodel & ~regTT.isConfPrev);

confparamsLT = nan([size(corr_mat_sess,1),size(regLT.confcoeff{idcLT},2)])
confparamsLT(1:90,:)  =  regLT.confcoeff{idcLT};
confparamsTT = nan([size(corr_mat_sess,1),size(regTT.confcoeff{idcTT},2)])
confparamsTT(1:90,:) =  regTT.confcoeff{idcTT};


%% load behavioral results



% calib_meanLT = (nanmean(conf_all,2)-nanmean(corr_all,2))*100;
conf_meanLT = nanmean(nanmean(nanmean(conf_mat_sess,2),3),4).*100;
acc_meanLT = nanmean(nanmean(nanmean(corr_mat_sess,2),3),4).*100;
calib_meanLT = conf_meanLT-acc_meanLT;

yLT = squeeze(nanmean(nanmean(conf_mat_sess,2),4));
valbias_meanLT = nanmean(yLT(:,[1,2])-yLT(:,[3,4]),2).*100;


% pref  = rlVars.pref;
% confpost = rlVars.confpost;
[prefMat_all,accMat_all,confMat_all,confMatPair_all,pref,accPost,confPost] = makePostMats(rlVars.ss,rlVars.aa,rlVars.confpost);

% accPost(accPost==0) = NaN;
conf_meanTT = nanmean(confPost,2).*100;
acc_meanTT = nanmean(accPost,2).*100;
calib_meanTT = conf_meanTT-acc_meanTT;
pair_mat = [6 5;2 1;8 7;4 3;];
yTT = [confPost(:,pair_mat(1,:)),confPost(:,pair_mat(2,:)),confPost(:,pair_mat(3,:)),confPost(:,pair_mat(4,:))];
valbias_meanTT = nanmean(yTT(:,[3,4,7,8])-yTT(:,[1,2,5,6]),2).*100;


acc_mean_confdataLT = acc_meanLT;
acc_mean_confdataLT(91:end,:) = NaN;

acc_mean_confdataTT = acc_meanTT;
acc_mean_confdataTT(91:end,:) = NaN;

allmeasLT = [acc_meanLT,conf_meanLT,valbias_meanLT,calib_meanLT];

allmeasTT = [acc_meanTT,conf_meanTT,valbias_meanTT,calib_meanTT];

params_learn = cell2mat(rlVars.params_allexp(:,nlearnmodel));

bias = (params_learn(:,2)-params_learn(:,3))./(params_learn(:,2)+params_learn(:,3));
bias2 = (params_learn(:,2)-params_learn(:,3));

learnparams_bias = [params_learn,bias,bias2];


npreg = load('Results\NonParRegSims.mat');% load non-par reg
% npreg = load('Results\NonParRegSimsBroadParams.mat');% load non-par reg


modelcolor = [.5,.5,.5;
    102,194,165;
    252,141,98;
    141,160,203]/255;

paramtypesnames = {'RL','conf'};
measnames = {'Accuracy','Confidence','Valence bias','Overconfidence'};
learnparamnames = {'$log(\beta)$','$\alpha_{CON}$','$\alpha_{DIS}$','$\alpha_V$','w','$\frac{\alpha_{CON}-\alpha_{DIS}}{\alpha_{CON}+\alpha_{DIS}}$','$\alpha_{CON}-\alpha_{DIS}$'};
confparamnames = {'Intercept','$\beta_{DIFF}$','$\beta_{BIAS}$','$\beta_{y(t-1)}$'};
for itask = 1:2
    switch itask
        case 1
            allmeas = allmeasLT;
            acc_confdata = acc_mean_confdataLT;
            myT = 'Learning';
            confparams = confparamsLT;
        case 2
            allmeas = allmeasTT;
            acc_confdata = acc_mean_confdataTT;

            myT = 'Transfer';
            confparams = confparamsTT;
    end
    for ipartype = 1:2
        switch ipartype
            case 1
                params = learnparams_bias;
                params(:,1) = log(params(:,1));
                paramnames = [learnparamnames];
            case 2
                params = confparams;
                paramnames = confparamnames;
        end
        
        R2 = cell(3,1);
        figure()
        [S,ax] = plotmatrix(params,allmeas)
        hold on
        title(myT)
        rData = nan(size(params,2),4);
        pData = nan(size(params,2),4);
        rTable = cell([size(params,2),4]);
        
        for ipar = 1:size(params,2)
            X = params(:,ipar);
            if ipartype ==1
                if ipar ==1
                    dumX = [0:0.1:max([X])]';
                elseif ipar==6
                    dumX = [-0.25:0.1:1]';
                elseif ipar>1
                    dumX = [0:0.1:1]';
                end
            elseif ipartype == 2
                dumX = [(min(X)-0.05*min(X)):0.1:(max(X)+0.05*max(X))]';
            end
            
            
            %%% accuracy for confidence data only
            Y = acc_confdata;
            [rData_acc_confexp(ipar),pData_acc_confexp(ipar)] = corr(X,Y,'type','Spearman','rows','complete');
            rDataAll_acc_confexp{itask,ipartype} = rData_acc_confexp;
            pDataAll_acc_confexp{itask,ipartype} = pData_acc_confexp;

            %%% all measures
            for imeas = 1:4
                set(S(imeas,ipar),'Marker','o','MarkerSize',4,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[0,0,0])
                Y = allmeas(:,imeas);
                [rData(ipar,imeas),pData(ipar,imeas)] = corr(X,Y,'type','Spearman','rows','complete');
                rDataAll{itask,ipartype} = rData;
                pDataAll{itask,ipartype} = pData;
                
                if ipartype == 1
                    if ipar ==1
                        xlim(ax(imeas,ipar),[min(dumX),max(dumX)])
                        xticks(ax(imeas,ipar),[0,1,2,3])
                    elseif ipar <6
                        if ipar == 3
                            xlim(ax(imeas,ipar),[0,.4])
                            xticks(ax(imeas,ipar),[0,0.2])
                        else
                            xlim(ax(imeas,ipar),[0,1])
                            xticks(ax(imeas,ipar),[0.2,0.5,0.8])
                        end
                    else
                        xlim(ax(imeas,ipar),[-0.25,1])
                        xticks(ax(imeas,ipar),[0,0.5])
                    end
                else
                    xlim(ax(imeas,ipar),[min(dumX),max(dumX)])
                    xlim(ax(imeas,ipar),[nanmean(X)-1.96*nanstd(X),nanmean(X)+1.96*nanstd(X)])
                end
                
                %xlim(ax(imeas,ipar),[min(dumX),max(dumX)])
                ylim(ax(imeas,ipar),[min(Y),max(Y)])
                
                %%%%%%%%% overlap non-parametric regressions %%%%%%
                
                for imodel = 2:4
                    if ~isempty(npreg.regAllX{imodel,itask,ipartype})
                        nprF = npreg.regAllF{imodel,itask,ipartype}{ipar,imeas};
                        nprX = npreg.regAllX{imodel,itask,ipartype}{ipar};
                        nprCI95 = npreg.regAllCI95{imodel,itask,ipartype}{ipar,imeas};
                        
                        l = line(ax(imeas,ipar),nprX,nprF,'Color',modelcolor(imodel,:),'LineWidth',3);
                        uistack(l,'bottom')
                        l = line(ax(imeas,ipar),nprX,nprCI95,'Color',modelcolor(imodel,:),'LineStyle','--');
                        uistack(l,'bottom')
                        
                        l =patch(ax(imeas,ipar),[nprX';flipud(nprX')],[nprCI95(2,:)';flipud(nprCI95(1,:)')],modelcolor(imodel,:),'LineStyle','None','FaceAlpha',0.2);
                        uistack(l,'bottom')
                        
                        
                        %goodness of fit
                        resid = nan(size(Y));
                        for iY = 1:numel(Y)
                            [~,idc] = min(abs(nprX-X(iY)));
                            resid(iY) = Y(iY)- nprF(idc);
                        end
                        ssRes = sum(resid.^2);
                        ssReg = sum((nprF-mean(Y)).^2);
                        varExp =  (ssReg)/((ssRes+ssReg));
                        R2{imodel}(ipar,imeas) = varExp;
                        text(ax(imeas,ipar),0.25,0.75,sprintf('\\rho = %0.2f',rData(ipar,imeas)),'Units','normalized','FontWeight','bold','FontSize',10,'BackgroundColor',[1,1,1,0.5])
                    end
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if ipar ==1
                    ylabel(ax(imeas,1),measnames{imeas},'Interpreter','latex')
                end
                xlabel(ax(4,ipar),paramnames{ipar},'Interpreter','latex')
                
                rTable{ipar,imeas} =['\rho$ = ',sprintf('%0.2f',rData(ipar,imeas)),...
                    ', p =',sprintf('%0.2f',pData(ipar,imeas))];
            end
        end
        
        
        
        dum = cell2table(rTable,'VariableNames',...
            {'Accuracy','Confidence','ValenceBias','Overconfidence'},...
            'RowNames',paramnames(1:size(params,2)));
        writetable(dum,['Results\corrBehavParams',myT,num2str(ipartype),'.csv'],...
            'WriteRowNames',true)
        
        for imodel = 1:3
            if ~isempty(R2{imodel})
                dum = array2table(round(R2{imodel},2),'VariableNames',...
                    {'Accuracy','Confidence','ValenceBias','Overconfidence'},...
                    'RowNames',paramnames(1:size(params,2)));
                writetable(dum,['Results\corrBehavParamsR2_MODEL',num2str(imodel),myT,num2str(ipartype),'.csv'],...
                    'WriteRowNames',true)
            end
        end
        
        switch ipartype
            case 1
                set(gcf,'Position',[0,0,800,500])
            case 2
                set(gcf,'Position',[0,0,500,500])
        end
        saveas(gcf,['Plots/corrBehavParams',paramtypesnames{ipartype},myT,'.svg'])
        
        
    end

end

save('Results/corrParamsBehav.mat','rDataAll','pDataAll','rDataAll_acc_confexp','pDataAll_acc_confexp');