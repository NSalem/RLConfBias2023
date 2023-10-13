
function genparams = getlearnparams(modelsinfo,n_subs,gen)
% get parameters for set of subjects and models for a single simulation
genparams = cell(n_subs,numel(modelsinfo));
for igenlearnmodel = 1:numel(modelsinfo)
    for iparam = 1:numel(modelsinfo{igenlearnmodel}.paramnames)
        thisParam = modelsinfo{igenlearnmodel}.paramnames{iparam};
        for isub = 1:n_subs
            genparams{isub,igenlearnmodel} = ...
                [genparams{isub,igenlearnmodel},gen.(thisParam)()];
        end
    end
end