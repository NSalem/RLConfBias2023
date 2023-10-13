addpath('helperfuncs')
load('Results\NonParRegSims.mat')


measnames = {'Accuracy','Confidence','Valence bias','Overconfidence'};
learnparamnames = {'$\beta$','$\alpha_{CON}$','$\alpha_{DIS}$','$\alpha_V$','w','$\frac{\alpha_{CON}-\alpha_{DIS}}{\alpha_{CON}+\alpha_{DIS}}$'};
modelnames = {'Q_C','\SigmaQ','V'};
confparamnames = {'$\beta_0$','$|\beta_{\Delta Q}|$','$\beta_{bias}$','$\beta_{y(t-1)}$'};
paramtypesnames = {'RL','conf'};
for imodel = 1:3
    for itask = 1:2
        if itask ==2 && imodel ==3
            break
        end
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
                    paramnames = confparamnames;
            end 

         for ipar =1:size(params,2)
             for imeas = 1:4
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
