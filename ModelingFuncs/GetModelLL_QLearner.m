function lik = GetModelLL_QLearner(params,learneroptions,s,a,r,c,aa,ss, usePriors, fitLearning,fitTransfer)
% Return the negative log likelihood or negative log posterior
% params: list of parameters
% learneroptions: structure with information about model 
% s: list of states in learning task
% a: list of actions in learning task 
% r: list of obtained outcomes in learning task 
% c: list of forgone outcomes in learning task
% ss: symbols presented each trial in transfer task
% aa: actions made in transfer task
% usePriors: whether to use priors to calculate the probability (meaning the function output will be log-posterior probability)
% fitLearning: whether to fit data on learning task 
% fitTransfer: whether to fit data on transfer task 

if nargin<11
    fitLearning = 1;
end
if nargin<10
    fitTransfer = 1;
end

%% Parameters
for iparam = 1:numel(learneroptions.paramnames)
    paramstruct.(learneroptions.paramnames{iparam}) = params(iparam);
end

for fn = fieldnames(learneroptions)'
    if ~strcmp(fn,'paramnames')
        paramstruct.(fn{1}) = learneroptions.(fn{1});
    end
end

l = qLearner(paramstruct);

lik     = 0;   % loglikelihood
for i = 1:length(a)
    if ~isnan(a(i))
        lik = lik+l.getActionLikelihood(s(i),a(i));
        
        COMP = (1-mod(s(i),2));
        PAR = round(mod(s(i),2));
        %%% learn
        if COMP
            l.learn(s(i),a(i),r(i),c(i));
        else
            l.learn(s(i),a(i),r(i));
        end
    end
end
Q = l.Q;

%% now post-test
Qf(1)=Q(end-4+1,2,end);
Qf(2)=Q(end-4+1,1,end);
Qf(3)=Q(end-4+2,2,end);
Qf(4)=Q(end-4+2,1,end);
Qf(5)=Q(end-4+3,2,end);
Qf(6)=Q(end-4+3,1,end);
Qf(7)=Q(end-4+4,2,end);
Qf(8)=Q(end-4+4,1,end);

PT_lik = 0;

for i = 1:length(aa)
    if (aa(i)) && (aa(i))~=1.5 % if a choice was performed in time
%         [action,p] = l.chooseActionPostLearning(Qf(ss(i,1)),Qf(ss(i,2))); % this might intruce small stochasticities due to rounding of p
                                                                          %it would be better to get the likelihood of the chosen option directly
% 
%         if action ~=aa(i)
%             p = 1-p;
%         end
%         PT_lik = PT_lik+log(p);

         PT_lik =PT_lik+l.getActionLikelihoodTransfer(Qf(ss(i,aa(i)))-Qf(ss(i,3-aa(i))));

    end
end
%%
lik = -lik.*fitLearning-PT_lik.*fitTransfer;

if usePriors
    lik = GetPosterior(l,lik,learneroptions.priorfuncs);
end
end

function [post]=GetPosterior(l,lik,priors)
%% Attribute probability to parameters
% used to calculate the LPP
% the priors we use are mostly taken from Daw et al. Neuron 2011

names = fieldnames(priors);
p = nan(1,numel(fieldnames(priors)));

for n  = 1:numel(names)
    iparam = names{n};
    ipriorfun = priors.(iparam);
    if ischar(class(ipriorfun))
        ipriorfun = str2func(ipriorfun);
    end
    p(n)= ipriorfun(l.(iparam));
end

p = -sum(p);

post = p + lik;
end
