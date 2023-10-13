%%%Make matrices for model identifiability and parameter
%recovery (Palminteri et al 2016,
%https://doi.org/10.1016/j.tics.2017.03.011)

clear
close all
clc
addpath('ModelingFuncs')
% load simulation data

simfiles = {'Results\SimsRLIdentRecov.mat'};
outfilenames = {'ResultsIdentFamilies.mat'};

for isimfile = 1
    fname = simfiles{isimfile};
    load(fname);
    
    if exist('isim','var') && isim~=nsims; %isim stores ongoing simulation if not finished
        nsims = isim-1;
    end
   
    nsubs = size(LAME,2);
    
    %%% identifiability for all family separations done in comparison on data %%%
    %1 SYM vs ASYM
    %2 ABS vs REL 
    %3 Contextualization types
    %4 ABS vs  REL vs CONF vs RELCONF
    familyGroups = {{[1:5],[6:10]},{[1,6],[2:5,7:10]},{[1,5],[2,6],[3,7],[4,8]},{1,[2:5],6,[7:10]}};
%     familyGroups = {{[1:5],[8:12]},{[1,8],[2:5,9:12]},{[2,9],[3,10],[4,11],[5,12]},{1,[2:5],8,[9:12]}};
    whichmodels = {[1:5,8:12],[1:5,8:12],[2:5,9:12],[1:5,8:12]};
    
    epAll = cell(size(familyGroups));
    bmAll = cell(size(familyGroups));
  
    for icomp = 1:numel(familyGroups);
        criterion = LAME;
        [epAll{icomp},bmAll{icomp}] = doModelIdent(LAME,whichmodels{icomp},familyGroups{icomp});
    end
 
    infostr = struct('nsim',nsims,'whichmodels',whichmodels); 
    save(['Results/',outfilenames{isimfile}],'bmAll','epAll','infostr')
end

