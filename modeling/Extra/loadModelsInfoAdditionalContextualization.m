modelsinfo{1} = struct('confirmatory',1,'contextual',1,'RuHatType','ROther','wcouType',2);
modelsinfo{1}.paramnames ={'beta','lr1', 'lr2','lr3','ROtherWeight'};
modelsinfo{1}.lb = [0,0,0,0,0];
modelsinfo{1}.ub = [Inf,1,1,1,1];
modelsinfo{1}.x0 = [1,.5,.5,.5,.5];

modelsinfo{2} = struct('confirmatory',1,'contextual',1,'RuHatType','ROther2');
modelsinfo{2}.paramnames ={'beta','lr1', 'lr2','lr3','ROtherWeight'};
modelsinfo{2}.lb = [0,0,0,0,0];
modelsinfo{2}.ub = [Inf,1,1,1,1];
modelsinfo{2}.x0 = [1,.5,.5,.5,.5];

modelsinfo{3} = struct('confirmatory',1,'contextual',1,'RuHatType','ROther');
modelsinfo{3}.paramnames ={'beta','lr1', 'lr2','lr3','ROtherWeight','fi','tau'};
modelsinfo{3}.lb = [0,0,0,0,0,-Inf,0];
modelsinfo{3}.ub = [Inf,1,1,1,1,Inf,1];
modelsinfo{3}.x0 = [1,.5,.5,.5,.5,0,.5];


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