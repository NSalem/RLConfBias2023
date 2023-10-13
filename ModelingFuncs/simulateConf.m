function [confL,confT] = simulateConf(rlvars,cfgL,cfgT)
%simulate confidence for a single participant
%inputs
% rlvars: struct/table with variables for confidence
% fields should correspond with names in cfg.varnames
% those for LT have dimensions sessions x trials
% for TT, vectors of trials
% cfg: struct
%       varnames: variables to include in confidence
%       coeff: coefficients for each variable
%       noise: sd for zero-mean additive gaussian noise
%outputs
%confLT: confidence in learning task (sessions x trials)
%confTT: confidence in transfer task (trials)


%confidence in learning task
confL = zeros(size(rlvars.dQ));
for isess = 1:size(rlvars.dQ,1)
    for itrial = 1:size(rlvars.dQ,2)
%         confLT(isess,itrial) = 0; 
        for ivar = 1:numel(cfgL.varnames)
            if strcmp(cfgL.varnames{ivar},'intercept')
                confL(isess,itrial) =  confL(isess,itrial)+cfgL.coeff(ivar);       
            elseif strcmp(cfgL.varnames{ivar},'confPrev') && itrial>1
                confL(isess,itrial) = confL(isess,itrial) + cfgL.coeff(ivar).*confL(isess,itrial-1);       
            elseif ~strcmp(cfgL.varnames{ivar},'confPrev')
                confL(isess,itrial) = confL(isess,itrial)+ cfgL.coeff(ivar).*rlvars.(cfgL.varnames{ivar})(isess,itrial); 
            end
        end
        
        confL(isess,itrial) = confL(isess,itrial) + randn*cfgL.noise;
        confL(isess,itrial) = 1./(1+exp(-confL(isess,itrial)));
%         confLT(isess,itrial) = 0.5+ 0.5/(1+exp(-confLT(isess,itrial)));
    end 
end

%confidence in transfer task
confT = zeros(size(rlvars.dQ_post));
for itrial = 1:size(rlvars.dQ_post,1)
%     confT(itrial) = 0;
    for ivar = 1:numel(cfgT.varnames)
        if strcmp(cfgT.varnames{ivar},'intercept')
            confT(itrial) = confT(itrial)+ cfgT.coeff(ivar);
        elseif strcmp(cfgT.varnames{ivar},'confPrev')&& itrial>1 
            confT(itrial) = confT(itrial) + cfgT.coeff(ivar)*confT(itrial-1);
        elseif ~strcmp(cfgT.varnames{ivar},'confPrev') 
            confT(itrial) = confT(itrial)+ cfgT.coeff(ivar).*rlvars.([cfgT.varnames{ivar},'_post'])(itrial);
        end 
    end
    
    confT(itrial) = confT(itrial) + randn*cfgT.noise;
    confT(itrial) = 1./(1+exp(-confT(itrial))); 
end
