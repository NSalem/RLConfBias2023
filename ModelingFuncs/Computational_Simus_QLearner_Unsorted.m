function [Q,V,pc,a,r,c,chpost, dQ_post,V_post,Q_post,Qf,Vf] = Computational_Simus_QLearner_Unsorted(params,s,ss,OUT)
%simulate choice and hidden variables given parameters,states, and potential outcomes
%INPUTS
%params: structure to be pased to qLearner
%s:  vector of context in learning task per trial (size trials)
%ss: matrix of symbols per trial in transfer task (2 x transfer trials)
%OUT: matrix of potential outcomes (options x trials)
%
%OUTPUTS
%Q: action values across learning task
%V: context values across learning task 
%pc: probability of correct choice across learning task
%r: obtained outcomes across learning task
%c: forgone outcomes across learning task
%chpost: choice over transfer task
%dQ_post: difference in presented values in transfer task
%V_post: average contextual value in transfer task
%Qf: Q-values for each symbol in transfer task (final values from learning)
%Vf: V-values for each symbol in transfer task  (final values from learning)

l = qLearner(params);
l.Q = zeros(numel(unique(s)),2);
l.V = zeros(numel(unique(s)));
%% Hidden variables
Q       = zeros(2,numel(s));        % Initial option values (all Models) as a function of condition ("s")
V       = zeros(1,numel(s));        % Initial Context values (Models 3) as a function of condition ("s")

c       = NaN(length(s(:)),1);
a       = NaN(length(s(:)),1);
pc      = NaN(length(s(:)),1);
r       = NaN(length(s(:)),1);
chpost  = NaN(length(ss),1);
dQ_post = NaN(length(ss),1);
Q_post = NaN(2,length(ss));
V_post = NaN(1,length(ss));

%%%% loop through trials %%%%
for i = 1:length(a)
    Q(:,i) = l.Q(s(i),:);
    V(i) = l.V(s(i));
    %         trialc(s(i)) =  trialc(s(i)) +1;
    %%% choose action
    [a(i),this_pc]=l.chooseAction(s(i));
    if a(i)~=2
        this_pc = 1-this_pc; %pc is probability of CORRECT choice
    end
    pc(i) = this_pc;
    
    %%% determine factual and counterfactual outcome
    r(i) = OUT(a(i),i);
    %         r(i) = OUT(a(i),i);
    c(i) = OUT(3-a(i),i);
    %          c(i) = OUT(3-a(i),i);
    COMP = (1-mod(s(i),2));
    PAR = round(mod(s(i),2));
    %%% learn
    if COMP
        l.learn(s(i),a(i),r(i),c(i));
    else
        l.learn(s(i),a(i),r(i));
    end
    
end

%% now post-test
nstates = numel(unique(s));

Qf(1)=l.Q(1,2); %PG75
Qf(2)=l.Q(1,1); %PG25
Qf(3)=l.Q(2,2); %CG75
Qf(4)=l.Q(2,1); %CG25
Qf(5)=l.Q(3,2); %PL25
Qf(6)=l.Q(3,1); %PL75
Qf(7)=l.Q(4,2); %CL25
Qf(8)=l.Q(4,1); %CL75


Vf(1)=l.V(1,1); %PG75
Vf(2)=l.V(1,1); %PG25
Vf(3)=l.V(2,1); %CG75
Vf(4)=l.V(2,1); %CG25
Vf(5)=l.V(3,1); %PL25
Vf(6)=l.V(4,1); %PL75
Vf(7)=l.V(4,1); %CL25
Vf(8)=l.V(4,1); %CL75


for k = 1: length(ss)
    dQ_post(k) = Qf(ss(k,2)) - Qf(ss(k,1));
    V_post(k) = (V(ss(k,2)) + Qf(ss(k,1)))/2;
    Q_post(:,k) = Qf(ss(k,:));
    
    chpost(k,1) = l.chooseActionPostLearning(Qf(ss(k,1)),Qf(ss(k,2)));
end
end
