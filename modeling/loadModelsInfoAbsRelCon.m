% This script creates a cell array of structs, each defining the
% specifications of a different model (e.g. whether to use contextual and
% confirmatory learning, which parameters does the model have, what are the
% bounds and prior distributions for each parameter)

modelsinfo{1}.counterfactual = 1; %ABS
modelsinfo{1}.paramnames ={'beta','lr1', 'lr2'};
modelsinfo{1}.lb = [0,0,0];
modelsinfo{1}.ub = [Inf,1,1];
modelsinfo{1}.x0 = [1,.5,.5];

modelsinfo{2} = struct('contextual',1,'counterfactual',1,'RuHatType','None'); %RELnone
modelsinfo{2}.paramnames = {'beta','lr1','lr2', 'lr3'};
modelsinfo{2}.lb = [0.,0,0,0];
modelsinfo{2}.ub = [Inf,1,1,1];
modelsinfo{2}.x0 = [1,.5,.5,.5];

modelsinfo{3} = struct('contextual',1,'counterfactual',1,'RuHatType','Qu');
modelsinfo{3}.paramnames = {'beta','lr1','lr2', 'lr3'};
modelsinfo{3}.lb = [0.,0,0,0];
modelsinfo{3}.ub = [Inf,1,1,1];
modelsinfo{3}.x0 = [1,.5,.5,.5];

modelsinfo{4} = struct('contextual',1,'counterfactual',1,'RuHatType','ROther');
modelsinfo{4}.paramnames ={'beta','lr1', 'lr2','lr3','ROtherWeight'};
modelsinfo{4}.lb = [0,0,0,0,0];
modelsinfo{4}.ub = [Inf,1,1,1,1];
modelsinfo{4}.x0 = [1,.5,.5,.5,.5];

modelsinfo{5} = struct('contextual',1,'counterfactual',1,'RuHatType','RLast');
modelsinfo{5}.paramnames = {'beta','lr1','lr2', 'lr3'};
modelsinfo{5}.lb = [0.,0,0,0];
modelsinfo{5}.ub = [Inf,1,1,1];
modelsinfo{5}.x0 = [1,.5,.5,.5];

modelsinfo{6} = struct('contextual',1,'counterfactual',1,'RuHatType','V');
modelsinfo{6}.paramnames ={'beta','lr1', 'lr2','lr3'};
modelsinfo{6}.lb = [0,0,0,0,0];
modelsinfo{6}.ub = [Inf,1,1,1,1];
modelsinfo{6}.x0 = [1,.5,.5,.5,.5];
modelsinfo{6}.RuHatType = 'V';

modelsinfo{7} = struct('contextual',1,'counterfactual',1,'RuHatType','V+Qu');
modelsinfo{7}.paramnames ={'beta','lr1', 'lr2','lr3'};
modelsinfo{7}.lb = [0,0,0,0,0];
modelsinfo{7}.ub = [Inf,1,1,1,1];
modelsinfo{7}.x0 = [1,.5,.5,.5,.5];


%%  Confirmatory
modelsinfo{8}.confirmatory = 1;
modelsinfo{8}.paramnames ={'beta','lr1', 'lr2'};
modelsinfo{8}.lb = [0,0,0];
modelsinfo{8}.ub = [Inf,1,1];
modelsinfo{8}.x0 = [1,.5,.5];

%% Confirmatory-Contextual
modelsinfo{9} = struct('confirmatory',1,'contextual',1,'RuHatType','None');
modelsinfo{9}.paramnames ={'beta','lr1', 'lr2','lr3'};
modelsinfo{9}.lb = [0,0,0,0];
modelsinfo{9}.ub = [Inf,1,1,1];
modelsinfo{9}.x0 = [1,.5,.5,.5];

modelsinfo{10} = struct('confirmatory',1,'contextual',1,'RuHatType','Qu');
modelsinfo{10}.paramnames ={'beta','lr1', 'lr2','lr3'};
modelsinfo{10}.lb = [0,0,0,0];
modelsinfo{10}.ub = [Inf,1,1,1];
modelsinfo{10}.x0 = [1,.5,.5,.5];
 
modelsinfo{11} = struct('confirmatory',1,'contextual',1,'RuHatType','ROther');
modelsinfo{11}.paramnames ={'beta','lr1', 'lr2','lr3','ROtherWeight'};
modelsinfo{11}.lb = [0,0,0,0,0];
modelsinfo{11}.ub = [Inf,1,1,1,1];
modelsinfo{11}.x0 = [1,.5,.5,.5,.5];

modelsinfo{12} = struct('confirmatory',1,'contextual',1,'RuHatType','RLast');
modelsinfo{12}.paramnames ={'beta','lr1', 'lr2','lr3'};
modelsinfo{12}.lb = [0,0,0,0];
modelsinfo{12}.ub = [Inf,1,1,1];
modelsinfo{12}.x0 = [1,.5,.5,.5];

modelsinfo{13} = struct('confirmatory',1,'contextual',1,'RuHatType','V');
modelsinfo{13}.paramnames ={'beta','lr1', 'lr2','lr3'};
modelsinfo{13}.lb = [0,0,0,0,0];
modelsinfo{13}.ub = [Inf,1,1,1,1];
modelsinfo{13}.x0 = [1,.5,.5,.5,.5];

modelsinfo{14} = struct('confirmatory',1,'contextual',1,'RuHatType','V+Qu');
modelsinfo{14}.paramnames ={'beta','lr1', 'lr2','lr3'};
modelsinfo{14}.lb = [0,0,0,0,0];
modelsinfo{14}.ub = [Inf,1,1,1,1];
modelsinfo{14}.x0 = [1,.5,.5,.5,.5];

priorfuncs.beta = '@(x) log(gampdf(x,1.2,5))';
priorfuncs.betaTT = '@(x) log(gampdf(x,1.2,5))';
priorfuncs.lr1 = '@(x) log(betapdf(x,1.1,1.1))';
priorfuncs.lr2 = priorfuncs.lr1; 
priorfuncs.lr3 = priorfuncs.lr1;    
priorfuncs.wcou = '@(x) log(betapdf(x/2,1.1,1.1))';
priorfuncs.ROtherWeight = '@(x) log(betapdf(x,1.1,1.1))';
priorfuncs.tau = priorfuncs.lr1;
priorfuncs.omega = priorfuncs.lr1;
priorfuncs.fi = '@(x) log(normpdf(0,1))';

for imodel = 1:numel(modelsinfo)
    for iparam = 1:numel(modelsinfo{imodel}.paramnames)
        thisParamName = modelsinfo{imodel}.paramnames{iparam};
        modelsinfo{imodel}.priorfuncs.(thisParamName) = priorfuncs.(thisParamName);
    end
end
