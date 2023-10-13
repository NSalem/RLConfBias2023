% Compare linear and logistic models of confidence

%% learning task 
reg = load('Results\reg_conflogit_learning_dqabs.mat');

BICL = [reg.BICCONF_linear,reg.BICCONF];
options.families = {1:12,13:24};
options.DisplayWin = 0;
[post,out] =VBA_groupBMC(-BICL'/2,options)

idc = 8;
BICL = [reg.BICCONF_linear(:,idc),reg.BICCONF(:,idc)];
[out,post] =VBA_groupBMC(-BICL'/2)

%% transfer task 
reg = load('Results\reg_conflogit_posttest_dqabs.mat');

BICT = [reg.BICCONF_linear,reg.BICCONF];
[post,out] = VBA_groupBMC(-BICT'/2,options);

idc = 2;
BICT = [reg.BICCONF_linear(:,idc),reg.BICCONF(:,idc)];
[out,post] =VBA_groupBMC(-BICT'/2)

