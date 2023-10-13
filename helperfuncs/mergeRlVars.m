
%%%% make mat concatenating all behav measures and hidden variables.
%%%% for learning task average across sessions


%%%% concatenate LT
rlvars_c.Q(:,:,4,:,:,:) = NaN;
rlvars_all.Q = cat(2,[rlvars_c.Q,rlvars_nc.Q]);

rlvars_c.correct(:,4,:,:) = NaN;
rlvars_all.correct = cat(2,[rlvars_c.correct,rlvars_nc.correct]);

varnames = {'dQ','dQabs','Q_c','Q_uc','V','pc','pa','conf','confprev'};
for ivar=varnames
    rlvars_c.(ivar{1})(:,:,4,:,:) = NaN;
    rlvars_all.(ivar{1}) = cat(2,[rlvars_c.(ivar{1}),rlvars_nc.(ivar{1})]);
end   

%%%% concatenate TT
varnames = {'dQ_post','dQabs_post','Q_c_post','Q_uc_post',...
    'V_post','pc_post','confpost','QG_post'};
for ivar=varnames
    rlvars_c.(ivar{1})(:,:,4,:,:) = NaN;
    rlvars_all.(ivar{1}) = cat(2,[rlvars_c.(ivar{1}),rlvars_nc.(ivar{1})]);
end   
    