function [pxpMat,bestMat] = doModelIdent(criterion,whichmodels,families)
%%% make model identifiability matrix
%%% criterion: matrix of size simulation x subject x generative model x
%%% fitted model
%%% pxpMat: matrix of size generative model x fitted model x simulation
%%% containing the protected exceedance probabilities for each model
%%% fitted on data from each generative model, for each simulation
%%%bestMat: matrix of size generative models x fitted models
%%% containing the number of simulations where each fitted model had the
%%% highest pxp for each generative model
if nargin<2
    whichmodels = 1:size(criterion,3);
end

if nargin<3
    families = [];
end

nmodels = numel(whichmodels);
nsims = size(criterion,1);

if ~isempty(families)
    bm      = zeros(numel(families),numel(families),nsims);  % best model
    pxpMat  = zeros(numel(families),numel(families),nsims);  % Protected exceedance probability
else
    bm      = zeros(nmodels,nmodels,nsims);  % best model
    pxpMat  = zeros(nmodels,nmodels,nsims);  % Protected exceedance probability
end
if ~isempty(families)
    options.families = families;
end

options.DisplayWin = 0;
for isim = 1:nsims
    for igenmodel = 1:numel(whichmodels)
        genmodel = whichmodels(igenmodel);
        [~, out] = VBA_groupBMC(squeeze(-criterion(isim,:,genmodel,whichmodels))',options);
        
        if ~isempty(families) 
            
            [~,winningfam] = max(out.families.ep);
            for ifam = 1:numel(families)
                if any(families{ifam}==igenmodel)
                    pxpMat(ifam,:,isim) = out.families.ep;
                    bm(ifam,winningfam,isim) = 1;
                end
            end
            
            
        else
            [~,winningmodel] = max(out.pxp);
            pxpMat(igenmodel,:,isim) = out.pxp;
            bm(igenmodel,winningmodel,isim) = 1;
        end
    end
    
    bestMat = sum(bm,3);
end
