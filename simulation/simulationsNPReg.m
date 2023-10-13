% Analyze the effect of parameters (learning and confidence) on
% simulated behavior using kernel smoothing regression, and save the
% results to a file. 

%We use the results of simulations based on parameters sampled from our
%(fitted) distributions in human data. For each RL parameter and confidence
%coefficient, we obtain a smooth line of predicted values for each
%behavioral measure (accuracy, confidence, overconfidence, valencee bias)
%for both learning and transfer task.


addpath('helperfuncs')
load('Results\SimulationsSubjParams.mat')
% load('Results\SimulationsSubjParamsUnif.mat')


ss_all = squeeze(outRL.ss);
ss_all = permute(ss_all,[1,3,2]);
aa_all = squeeze(outRL.choicePost);

thispc = squeeze(nanmean(nanmean(nanmean(outRL.pc,3),2),4));
for imodel = 1:size(outConfLSorted,3)
    meanConfCond = squeeze(nanmean(nanmean(nanmean(outConfLSorted(:,:,imodel,:,:,:),1),4),6));
    confBiasLT{imodel} = (squeeze(mean(meanConfCond(:,[1,2])-meanConfCond(:,[3,4]),2)))*100;
    overconfLT{imodel} = mean(meanConfCond,2)*100-thispc.*100;
    confLT{imodel} = squeeze(mean(meanConfCond,2)*100);
    
end

newIdc = [6,5,2,1,8,7,4,3]; %rearange order of post-test stimuli
accTT = squeeze(nanmean(accPost(:,:,newIdc))).*100;
confPost = squeeze(nanmean(confPost(:,:,:,newIdc)))*100;

for imodel= 1:size(outConfT,3)
    confTT = squeeze(confPost(imodel,:,:));
    confTTAll{imodel} = confTT;
    overconfTT{imodel} = nanmean(confTT-accTT,2);
    confBiasTT{imodel} = nanmean(confTT(:,[3,4,7,8])-confTT(:,[1,2,5,6]),2);
end

learnbias = (learnparams(:,2)-learnparams(:,3))./((learnparams(:,2)+learnparams(:,3)));
learnbias2 = (learnparams(:,2)-learnparams(:,3));
learnpars_bias = [learnparams,learnbias,learnbias2];

for imodel = 1:4
    allparsLT{imodel} = [learnparams,confparamsLT{imodel}];
    allmeasLT{imodel} = [thispc.*100, confLT{imodel},confBiasLT{imodel},overconfLT{imodel}];
    allparsTT{imodel} = [learnparams,confparamsTT{imodel}];
    allmeasTT{imodel}  = [nanmean(accTT,2),nanmean(confTTAll{imodel},2),confBiasTT{imodel},overconfTT{imodel}];
end

for itask = 1:2
    switch itask
        case 1
            allmeas = allmeasLT{1};
            confparams = confparamsLT{1};
            myT = 'Learning';
        case 2
            myT = 'Transfer';
            allmeas = allmeasTT{1};
            confparams = confparamsTT{1};
    end
end

measnames = {'Accuracy','Confidence','Valence bias','Overconfidence'};
learnparamnames = {'$\beta$','$\alpha_{CON}$','$\alpha_{DIS}$','$\alpha_V$','w','$\frac{\alpha_{CON}-\alpha_{DIS}}{\alpha_{CON}+\alpha_{DIS}}$','$\alpha_{CON}-\alpha_{DIS}$'};
modelnames = {'unbiased','Q_C','\SigmaQ','V'};
confparamnames = {'$\beta_0$','$|\beta_{\Delta Q}|$','$\beta_{bias}$','$\beta_{y(t-1)}$'};
paramtypesnames = {'RL','conf'};
for imodel = 1:4
    for itask = 1:2
        switch itask
            case 1
                confpars = confparamsLT{imodel};
                allmeas = allmeasLT{imodel};
                taskname = 'Learning';
            case 2
                confpars = confparamsTT{imodel};
                allmeas = allmeasTT{imodel};
                taskname = 'Transfer';
        end
        myT = [taskname,', ' ,modelnames{imodel}];
        
        %%% scatterplot matrix of learning params vs behavior
        for ipartype = 1:2            
            switch ipartype
                case 1
                    params = learnpars_bias;
                    params(:,1) = log(params(:,1));
                    paramnames = learnparamnames;
                case 2
                    params = confpars;
                    %                     params = 1./(1+exp(-params));
                    paramnames = confparamnames;
            end
            
            rTableModel = cell([size(params,2),4]);
            
            for ipar = 1:size(params,2)
                X = params(:,ipar);
                for imeas = 1:4
                    Y = allmeas(:,imeas);
                    [rho{imodel,itask,ipartype}(ipar,imeas),p{imodel,itask,ipartype}(ipar,imeas)] = corr(X,Y,'type','Spearman');
                   
                    %%% parametric regression
                    dum = linspace(min(X),max(X),1000)';
                    if strcmp(paramnames{ipar},'$\alpha_V$')
                        thislm = fitlm(log(X),Y,'RobustOpts','on');
                        [y,yci] = thislm.predict(log(dum));
                        
                    else
                        thislm = fitlm(X,Y,'RobustOpts','on');
                        [y,yci] = thislm.predict(dum);
                    end
                    
                    regLin{imodel,itask,ipartype}{ipar,imeas} = thislm;
                                        
                    %%% non parametric regression (gaussian kernel
                    %%% smoothing)
                    %%%% choosing an h with cross-validation
                    hVals = linspace(range(X)/50,range(X),50);
                    CVerr = nan(size(hVals));
                    for ih=1:numel(hVals)
                        %%% cross-validation folds
                        CVerr(ih) = cvKfoldKsr(X,Y,5,hVals(ih),100);
                    end
                    [~,idmin]= min(CVerr);
                    
                    npreg = ksr(X,Y,hVals(idmin));
                    
                    %%% make ci with bootsrapping
                    %%% do 100 resamples of 100 datapoints and make the
                    %%% regression, then calculate ci95 point by point
                    
                    nResamples = 100;
                    nPointsSample = 100;
                    fStar = nan(nResamples,numel(npreg.f));
                    for iresample = 1:nResamples
                        idc = randi([1,numel(X)],1,nPointsSample);
                        xStar = X(idc);
                        yStar = Y(idc);
                        regBootstrap = ksr(xStar,yStar,npreg.h,npreg.N,[min(npreg.x),max(npreg.x)]);
                        fStar(iresample,:) = regBootstrap.f;
                        clear regBootstrap;
                    end
                    
                    npreg.ci95 =  [npreg.f+1.96.*std(fStar);npreg.f-1.96.*std(fStar)];
                    regAllF{imodel,itask,ipartype}{ipar,imeas} = npreg.f;
                    %                     regAllSD{imodel,itask,ipartype}{ipar,imeas} = npreg.sd;
                    regAllX{imodel,itask,ipartype}{ipar} = npreg.x;
                    regAllCI95{imodel,itask,ipartype}{ipar,imeas} = npreg.ci95;
                    
                    %goodness of fit
                    resid = nan(size(Y));
                    for iY = 1:numel(Y)
                        [~,idc] = min(abs(npreg.x-X(iY)));
                        resid(iY) = Y(iY)- npreg.f(idc);
                    end
                    ssRes = sum(resid.^2);
                    ssReg = sum((npreg.f-mean(Y)).^2);
                    varExp =  (ssReg)/((ssRes+ssReg));
                    
                    rTableModel{ipar,imeas}=['\rho$ = ',sprintf('%0.2f',rho{imodel,itask,ipartype}(ipar,imeas)),...
                        ', p =',sprintf('%0.2f',p{imodel,itask,ipartype}(ipar,imeas))];
                end
            end
            
            
            dum = cell2table(rTableModel,'VariableNames',...
                {'Accuracy','Confidence','ValenceBias','Overconfidence'},...
                'RowNames',paramnames(1:size(params,2)));
            writetable(dum,['Results\corrBehavParamsMODEL',num2str(imodel),taskname,num2str(ipartype),'.csv'],...
                'WriteRowNames',true)
        end
    end
end
save('Results\NonParRegSims.mat','regAllX','regAllF','regAllCI95','rho', 'p','regLin');


rmpath('helperfuncs')
rmpath('ModelingFuncs')