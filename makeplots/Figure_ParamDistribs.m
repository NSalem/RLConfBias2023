%Plot distributions (pirate plots) of RL parameters and confidence coefficients as fitted on participants' behavior

addpath('helperfuncs')

%% RL parameters for best fitting model  (RELASYMOther)

rl = load('Results\RLvars_all.mat');
imodel = 11; %RELASYMOther
 
params = [cell2mat(rl.params_allexp(:,imodel))];

figure();
col = repmat(.5,size(params,2),3);
subplot(1,4,1)
pirateplot(params(:,1)',col(1,:),0,30,12,'','','','')
xticklabels('\beta');

subplot(1,4,[2:4])
pirateplot(params(:,2:end)',col,0,1,12,'','','','')
xticklabels({'\alpha_{CON}','\alpha_{DIS}','\alpha_V','w'})

%prior distributions (check formulas match modelsinfo{imodel}.priors;
for ipar = 1:size(params,2)
    if ipar == 1
        p_thispar = gamfit(params(:,ipar));
        
    else
        p_thispar = betafit(params(:,ipar));
        
    end
    p{ipar} = p_thispar;
end

x = linspace(0,100,1000);
x2 = linspace(0,1,1000);
priorBeta = gampdf(x,1.25,5);
priorLR = betapdf(x2,1.1,1.1);
priorW = betapdf(x2,1.1,1.1);
priorf=cell(5,1);
priorf{1} = priorBeta;

distf = cell(5,1);
for ipar = 1:5
    if ipar ==1
        distf{ipar} = gampdf(x,p{ipar}(1),p{ipar}(2));
    else
        distf{ipar} = betapdf(x2,p{ipar}(1),p{ipar}(2));
    end
end

for k = 2:4
    priorf{k} = priorLR;
end
priorf{5} = priorW;

for iparam = 1:numel(rl.modelsinfo{imodel}.paramnames)
    width(iparam) = 0.75/2/max(priorf{iparam}); %to match pirateplot
    widthdist(iparam) = 0.75/2/prctile(distf{iparam},99); %to match pirateplot
end

priorcol = [0.2,0.4,1];
distcol = [0.7,0,0];
subplot(1,4,1)
plot(priorBeta*width(1)+1,x,'LineStyle','--','Color',priorcol);
plot(-priorBeta*width(1)+1,x,'LineStyle','--','Color',priorcol);

plot(distf{1}*widthdist(1)+1,x,'LineStyle','--','Color',distcol);
plot(-distf{1}*widthdist(1)+1,x,'LineStyle','--','Color',distcol);


subplot(1,4,[2:4])
hold on
% pirateplot(params(:,2:5)',col(2:5,:),0,1,12,'','','','')

for iparam = 2:numel(rl.modelsinfo{imodel}.paramnames)
    plot(priorf{iparam}*width(iparam)+iparam-1,x2,'LineStyle','--','Color',priorcol);
    plot(-priorf{iparam}*width(iparam)+iparam-1,x2,'LineStyle','--','Color',priorcol);
    plot(distf{iparam}*widthdist(iparam)+iparam-1,x2,'LineStyle','--','Color',distcol);
    plot(-distf{iparam}*widthdist(iparam)+iparam-1,x2,'LineStyle','--','Color',distcol);

end

saveas(gcf,['Plots/paramsRL.svg']);

% xticklabels({'\alpha_{CON}','\alpha_{DIS}','\alpha_V','w'})

clear;


%% Confidence parameters (coefficients) for all bias models (Qc,sigmaQ,V) and both tasks (learning and transfer)
regLT = load('Results\reg_conflogit_learning_dqabs.mat'); %learning task
idcLT = [8,9,10];%find(regLT.whichConfModel == 2 & regLT.whichLearnModel == 11 & regLT.isConfPrev);

regTT = load('Results\reg_conflogit_posttest_dqabs.mat'); %transfer task
idcTT = [2,3,4];%find(regTT.whichConfModel == 2 & regTT.whichLearnModel == 11 & ~regTT.isConfPrev);

labelsLT = {'B_{0}','B_{|\DeltaQ|}','B_{BIAS}','B_{Conf(t-1)}'};
labelsTT = {'B_{0}','B_{|\DeltaQ|}','B_{BIAS}'} ;
biasNames = {'B_{Q_C}','B_{\SigmaQ}','B_V'};
modelnames = {'Qc','simgaQ','V'}
% coeffLT(:,4) = coeffLT(:,4);

% coeffTT(:,4) = coeffTT(:,4)*100;

col = repmat(.5,4,3);
priorcol = [0.2,0.4,1];
distcol = [0.7,0,0];

for imodel = 1:3 %loop over confidence models (Qc bias, sigmaQ bias, V bias)
    coeffLT = regLT.confcoeff{idcLT(imodel)};
    coeffTT = regTT.confcoeff{idcTT(imodel)};
    labelsLT{3} = biasNames{imodel};
    labelsTT{3} = biasNames{imodel};
    
    figure()
    subplot(1,2,1)
    hold on;
    plot([0,5],[0,0],'k--')
    pirateplot(coeffLT',col,-3,8,12,'Learning Task','','');
 
    xticklabels(labelsLT);
    xtickangle(90);
    for iparam = 1:size(coeffLT,2)
        
        idc = all(abs(zscore(coeffLT))<=1.96,2);
        pm = mean(coeffLT(idc,iparam));
        ps = std(coeffLT(idc,iparam));
        
        plot(normpdf(-7:0.2:7,pm,ps)*0.75+iparam,-7:0.2:7,'LineStyle','--','Color',distcol);
        plot(-normpdf(-7:0.2:7,pm,ps)*0.75+iparam,-7:0.2:7,'LineStyle','--','Color',distcol);
    end
    
    xtickangle(90);
    subplot(1,2,2)
    hold on
    plot([0,5],[0,0],'k--')
    pirateplot(coeffTT',col,-3,8,12,'Transfer Task','','');
    xticklabels(labelsTT);
    xtickangle(90);
    for iparam = 1:size(coeffTT,2)
        
        idc = all(abs(zscore(coeffTT))<=1.96,2);
        pm = mean(coeffTT(idc,iparam));
        ps = std(coeffTT(idc,iparam));
        
        plot(normpdf(-7:0.2:7,pm,ps)*0.75+iparam,-7:0.2:7,'LineStyle','--','Color',distcol);
        plot(-normpdf(-7:0.2:7,pm,ps)*0.75+iparam,-7:0.2:7,'LineStyle','--','Color',distcol);
    end
    saveas(gcf,['Plots/confParams_model',modelnames{imodel},'.svg']);
end

rmpath('helperfuncs')
