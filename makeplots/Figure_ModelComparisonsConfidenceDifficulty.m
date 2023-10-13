
%%% make comparison between confidence models for different types of
%%% difficulty variables


addpath('helperfuncs')

files = {'Results\reg_conflogit_learning_dqabs','Results\reg_conflogit_learning','Results\reg_conflogit_learning_dqgoodbad','Results\reg_conflogit_learning_pc'}
% files = {'Results\reg_conflogit_norandom_learning_dqabs','Results\reg_conflogit_norandom_learning','Results\reg_conflogit_norandom_learning_dqgoodbad','Results\reg_conflogit_norandom_learning_pc'}

%% compare bias models for each difficulty model, for Learning, Transfer, and both tasks
figure()
BICall = [];
LLall= [];
nparsall = [];
for ip = 1:4
    load(files{ip})
    
    subplot(3,4,ip)
  
    [postLT(ip),outLT(ip)] = BMC_plot(BICCONF(:,[7,8,10]));
    
    xticklabels({'None','Q_C','V'})
    xlabel('bias')
    
    BICall = [BICall,BICCONF(:,[7,8,10])];
    LLall  =[LLall,LLCONF(:,[7,8,10])];
    nparsall =[nparsall,npars(:,[7,8,10])];
    clear BICCONF LLCONF npars
end


files = {'Results\reg_conflogit_posttest_dqabs','Results\reg_conflogit_posttest',...
    'Results\reg_conflogit_posttest_dqgoodbad','Results\reg_conflogit_posttest_pc'};

% 
% files = {'Results\reg_conflogit_norandom_posttest_dqabs','Results\reg_conflogit_norandom_posttest',...
%     'Results\reg_conflogit_norandom_posttest_dqgoodbad','Results\reg_conflogit_norandom_posttest_pc'};


BICallTT = [];
LLallTT = [];
nparsallTT = [];
for ip = 1:4
    load(files{ip})
    
    subplot(3,4,ip+4)
  
     [postTT(ip),outTT(ip)] = BMC_plot(BICCONF(:,[7,8,10]));
    
    xticklabels({'None','Q_C','V'})
    xlabel('bias')
    
    BICallTT = [BICallTT,BICCONF(:,[7,8,10])];
    LLallTT = [LLallTT,LLCONF(:,[7,8,10])];
    nparsallTT =[nparsallTT,npars(:,[7,8,10])]; 
    
end

LLallJoint = LLall+LLallTT;

BICallJoint = (nparsall+nparsallTT)*log(112+288)-2*LLallJoint;

families =  {[1,2,3],[4,5,6],[7,8,9],[10,11,12]};
for ip = 1:4
    
    subplot(3,4,ip+8)
  
    [postJoint(ip),outJoint(ip)] = BMC_plot(BICallJoint(:,families{ip}));
    
    xticklabels({'None','Q_C','V'})
    xlabel('bias')
  
    
end

set(gcf,'Position',[0,0,1000,600])
saveas(gcf,'Plots/BMCConfBias_AllDifficultyFamilies_withSigmaQ.svg')



options.families = {[1,2,3],[4,5,6],[7,8,9],[10,11,12]};
labels = {'|\DeltaQ|','Q_C-Q_U','Q_{GOOD}-Q_{BAD}','p(Correct)'};

figure
subplot(1,3,1)
[post,outLT] = BMC_plot(BICall, options)
xticklabels(labels)
subplot(1,3,2)
[post,outTT] = BMC_plot(BICallTT, options)
xticklabels(labels)

subplot(1,3,3)
[post,outJoint] = BMC_plot(BICallJoint, options)
xticklabels(labels)

set(gcf, 'Position',[0,0,1500,400])
saveas(gcf,'Plots/BMCConfDifficulty_withSimgaQ.svg')


% options.families = {[1,2,3],[4,5,6],[7,8,9]};
% labels = {'|\DeltaQ|','Q_C-Q_U','Q_{GOOD}-Q_{BAD}'};
% figure
% subplot(1,3,1)
% [post,outLT] = BMC_plot(BICall(:,[1:9]), options)
% xticklabels(labels)
% subplot(1,3,2)
% [post,outTT] = BMC_plot(BICallTT(:,[1:9]), options)
% xticklabels(labels)
% 
% subplot(1,3,3)
% [post,outJoint] = BMC_plot(BICallJoint(:,[1:9]), options)
% xticklabels(labels)
% 
% set(gcf, 'Position',[0,0,1500,400])



options.families = {[1,2,3],[4,5,6]};
labels = {'|\DeltaQ|','Q_C-Q_U'};
figure()
subplot(1,3,1)
[post,outLT] = BMC_plot(BICall(:,[1:6]), options)
xticklabels(labels)
subplot(1,3,2)
[post,outTT] = BMC_plot(BICallTT(:,[1:6]), options)
xticklabels(labels)

subplot(1,3,3)
[post,outJoint] = BMC_plot(BICallJoint(:,[1:6]), options)
xticklabels(labels)

set(gcf, 'Position',[0,0,1500,400])
saveas(gcf,'Plots/BMCConfDifficulty_SignedVsUnsigned_withSigmaQ.svg')


%% same but extended w/ sigmaq where possible (all diff models except qc-qu)
clear

files = {'Results\reg_conflogit_learning_dqabs','Results\reg_conflogit_learning','Results\reg_conflogit_learning_dqgoodbad','Results\reg_conflogit_learning_pc'}

figure()
BICall = [];
LLall= [];
nparsall = [];
for ip = 1:4
    load(files{ip})
    
    subplot(3,4,ip)

    
    if ip ~= 2
        ind = 7:10;
        xl = {'None','Q_C','\Sigma_Q','V'};
    else 
        ind = [7,8,10];
        xl = {'None','Q_C','V'};
    end
    
    [postLT(ip),outLT(ip)] = BMC_plot(BICCONF(:,ind));
    
    xticklabels(xl)
    xlabel('bias')
    
    BICall = [BICall,BICCONF(:,ind)];
    LLall  =[LLall,LLCONF(:,ind)];
    nparsall =[nparsall,npars(:,ind)];
end

files = {'Results\reg_conflogit_posttest_dqabs','Results\reg_conflogit_posttest',...
    'Results\reg_conflogit_posttest_dqgoodbad','Results\reg_conflogit_posttest_pc'};

BICallTT = [];
LLallTT = [];
nparsallTT = [];
for ip = 1:4
    load(files{ip})
      
    if ip ~= 2
        ind = 7:10;
        xl = {'None','Q_C','\Sigma_Q','V'};
    else 
        ind = [7,8,10];
        xl = {'None','Q_C','V'};
    end
     
    subplot(3,4,ip+4)
  
     [postTT(ip),outTT(ip)] = BMC_plot(BICCONF(:,ind));
    
    xticklabels(xl)
    xlabel('bias')
    
    BICallTT = [BICallTT,BICCONF(:,ind)];
    LLallTT = [LLallTT,LLCONF(:,ind)];
    nparsallTT =[nparsallTT,npars(:,ind)];
    
end

LLallJoint = LLall+LLallTT;

BICallJoint = (nparsall+nparsallTT)*log(112+288)-2*LLallJoint;

families = {[1,2,3,4],[5,6,7],[8,9,10,11],[12,13,14,15]};
for ip = 1:4
    
    if ip ~= 2
        xl = {'None','Q_C','\Sigma_Q','V'};
    else 
        xl = {'None','Q_C','V'};
        end
    
    subplot(3,4,ip+8)
  
    [postJoint(ip),outJoint(ip)] = BMC_plot(BICallJoint(:,families{ip}));
    
    xticklabels(xl)
    xlabel('bias')
end



set(gcf,'Position',[0,0,1000,600])
saveas(gcf,'Plots/BMCConfBias_AllDifficultyFamilies_withoutSigmaQ.svg')


options.families = families;
labels = {'|\DeltaQ|','Q_C-Q_U','Q_{GOOD}-Q_{BAD}','p(Correct)'};

figure
subplot(1,3,1)
[post,outLT] = BMC_plot(BICall, options)
xticklabels(labels)
subplot(1,3,2)
[post,outTT] = BMC_plot(BICallTT, options)
xticklabels(labels)

subplot(1,3,3)
[post,outJoint] = BMC_plot(BICallJoint, options)
xticklabels(labels)


set(gcf, 'Position',[0,0,1500,400])
saveas(gcf,'Plots/BMCConfDifficulty_withoutSigmaQ.svg')
% 
% 
% options.families = {[1,2,3,4],[5,6,7],[8,9,10,11]};
% labels = {'|\DeltaQ|','Q_C-Q_U','Q_{GOOD}-Q_{BAD}'};
% figure
% subplot(1,3,1)
% [post,outLT] = BMC_plot(BICall(:,[1:11]), options)
% xticklabels(labels)
% subplot(1,3,2)
% [post,outTT] = BMC_plot(BICallTT(:,[1:11]), options)
% xticklabels(labels)
% 
% subplot(1,3,3)
% [post,outJoint] = BMC_plot(BICallJoint(:,[1:11]), options)
% xticklabels(labels)
% 
% set(gcf, 'Position',[0,0,1500,400])
% 
% 

options.families = {[1,2,3,4],[5,6,7]};
labels = {'|\DeltaQ|','Q_C-Q_U'};
figure
subplot(1,3,1)
[post,outLT] = BMC_plot(BICall(:,[1:7]), options)
xticklabels(labels)
subplot(1,3,2)
[post,outTT] = BMC_plot(BICallTT(:,[1:7]), options)
xticklabels(labels)

subplot(1,3,3)
[post,outJoint] = BMC_plot(BICallJoint(:,[1:7]), options)
xticklabels(labels)

set(gcf, 'Position',[0,0,1500,400])
saveas(gcf,'Plots/BMCConfDifficulty_SignedVsUnsigned_withoutSigmaQ.svg')


%% compare y~|dQ|+Qc with y~Qc+Qu
% load agnostic models
regAbsDQ_L = load('Results\reg_conflogit_learning_dqabs.mat');
regAbsDQ_T = load('Results\reg_conflogit_posttest_dqabs.mat');

regNoDQ_L = load('Results\reg_conflogit_learning_nodq.mat');
regNoDQ_T =  load('Results\reg_conflogit_posttest_nodq.mat');
LLNoDQJoint = regNoDQ_L.LLCONF + regNoDQ_T.LLCONF;
BICNoDQJoint = (regNoDQ_L.npars+regNoDQ_T.npars)*log(112+288)-2*LLNoDQJoint;

dum_L = [regAbsDQ_L.BICCONF(:,8),regNoDQ_L.BICCONF(:,9)];    
dum_T = [regAbsDQ_T.BICCONF(:,8),regNoDQ_T.BICCONF(:,9)];
dum_J = [BICallJoint(:,8),BICNoDQJoint(:,9)];

figure()
options = {};
subplot(1,3,1)
[post,outLT] = BMC_plot(dum_L, options);
xticklabels({'|dQ|+Qc','Qc+Qu'})
subplot(1,3,2)
[post,outLT] = BMC_plot(dum_T, options);
xticklabels({'|dQ|+Qc','Qc+Qu'})
subplot(1,3,3)
[post,outLT] = BMC_plot(dum_J, options);
xticklabels({'|dQ|+Qc','Qc+Qu'})
