function [Rest,R2est] = doParamRecov(genpar,recpar)
%%% genpar: matrix of generative parameters size subj x param x sim
%%% recpar: matrix of fitted parameters size subj x param x sim

nparams =size(genpar,2);
nsims = size(genpar,3);
Rest    = NaN(nparams,nparams,nsims);
R2est   = NaN(nparams,nparams,nsims);

for isim = 1:nsims
    % compute correlations between parameters used to simulate the data,
    % and recovered (i.e. estimated) parameters
    Rest(:,:,isim) = corr(squeeze(recpar(:,:,isim)));
    for ipar = 1:nparams
        Rest(ipar,ipar,isim) = corr(squeeze(recpar(:,ipar,isim)),squeeze(genpar(:,ipar,isim)));
    end
    R2est(:,:,isim) = Rest(:,:,isim).*Rest(:,:,isim);
end

end