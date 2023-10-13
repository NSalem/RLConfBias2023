function [confL,confT,confLsorted] = simulateConfPop(genparamsL,genparamsT,noiseL,noiseT,modelsinfo,modelsinfoT,rlVars)
% generate confidence for learning and transfer task, for each model
%genparamsL/genparamsT are cell arrays of lenght n models, containing matrix of size
%participants x parameter

for iconfmodel = 1:numel(modelsinfo)
    ilearnmodel = modelsinfo{iconfmodel}.ilearnmodel;
    
    for isub = 1:size(genparamsL{1},1)
        paramssubL = genparamsL{iconfmodel}(isub,:);
        paramssubT = genparamsT{iconfmodel}(isub,:);
     
        %make rlvarsSub
        rlVarsSub = struct('dQ', squeeze(rlVars.dQ(isub,ilearnmodel,:,:)),...
            'dQabs', squeeze(rlVars.dQabs(isub,ilearnmodel,:,:)),...
            'Qc', squeeze(rlVars.Qc(isub,ilearnmodel,:,:)),...
            'Qu',squeeze(rlVars.Qu(isub,ilearnmodel,:,:)),...
            'V', squeeze(rlVars.V(isub,ilearnmodel,:,:)),...
            'sigmaQ', squeeze(rlVars.sigmaQ(isub,ilearnmodel,:,:)),...
            'dQ_post', squeeze(rlVars.dQ_post(isub,ilearnmodel,:)),...
            'dQabs_post', squeeze(rlVars.dQabs_post(isub,ilearnmodel,:)),...
            'Qc_post', squeeze(rlVars.Qc_post(isub,ilearnmodel,:)),...
            'Qu_post', squeeze(rlVars.Qu_post(isub,ilearnmodel,:)),...
            'sigmaQ_post', squeeze(rlVars.sigmaQ_post(isub,ilearnmodel,:)),...
            'V_post', squeeze(rlVars.V_post(isub,ilearnmodel,:)));
        
        %make cfgL and cfgT
        cfgLsub = struct();
        cfgLsub.varnames = modelsinfo{iconfmodel}.varnames;
        cfgLsub.coeff = paramssubL;
        cfgLsub.noise = noiseL{ilearnmodel}(isub);
        cfgTsub = struct();
        cfgTsub.varnames = modelsinfoT{iconfmodel}.varnames;
        cfgTsub.coeff = paramssubT;
        cfgTsub.noise = noiseT{iconfmodel}(isub);
        
        [confLsub, confTsub]= simulateConf(rlVarsSub,cfgLsub,cfgTsub);
    
        %%%  store data %%%
        confL(isub,iconfmodel,:,:) = confLsub;
        confT(isub,iconfmodel,:,:) = confTsub;
 
        nsess = size(confLsub,1);
        s = squeeze(rlVars.s(isub,1,:,:));
        for isess = 1:nsess
            for icond = 1:4
                confLsorted(isub,iconfmodel,isess,icond,:) = confLsub(isess,s(isess,:)==icond);
            end
        end 
    end
end
end
