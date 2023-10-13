function [Q,V,pc,PE,ppc,pll,dQ_post,V_post,Q_post,Qf,Vf] = Computational_TimeSeries_QLearner_Unsorted(params,s,a,r,c,aa,ss)
%simulate hidden variables given actions and parameters
%Inputs:
%params: struct containing model specifications and parameters (as input for creating instance of qLearner class)
%s: list of states in learning task
%a: list of actions in learning task
%r: list of obtained outcomes in learning task
%c: list of forgone outcomes in learning task
%ss: list of symbols presented in transfer task
%aa: list of actions in transfer task
%Outputs:
%Q: Action values over learning task
%V: context values over learning task
%pc: probability of correct choice during learning task
%ppc: probability of chosing symbol on the right during transfer task
%pll: likelihood of choice in transfer task

%% start qlearner
l = qLearner(params);
l.Q = zeros(numel(unique(s)),2);
l.V = zeros(numel(unique(s)));
ntrials = numel(s);
Q       = zeros(2,ntrials);        % Initial option values (all Models) as a function of conditio ("s")
V       = zeros(1,ntrials);        % Initial Context values (Models 3) as a function of conditio ("s")
pc      = ones(1,ntrials)*.5;
dQ_post = NaN(length(ss),1);
V_post = NaN(length(ss),1);
Q_post = NaN(2,length(ss));

for i = 1:length(a)
    Q(:,i) = l.Q(s(i),:);   
    V(i) = l.V(s(i));        
    if ~isnan(a(i))
        %%% get probability of action being correct given current values 
        [dum_action,trl_pc] = l.chooseAction(s(i)); 
        if dum_action ~=2
            trl_pc = 1-trl_pc; 
        end
        pc(i) = trl_pc;
        COMP = (1-mod(s(i),2));
        PAR = round(mod(s(i),2));
        %%% learn
        if COMP
            l.learn(s(i),a(i),r(i),c(i));
        else
            l.learn(s(i),a(i),r(i));
        end  
        
        if i<length(a)
            PE(i+1) = r(i)-l.V(s(i)) - l.Q(s(i),a(i));
        end        
    end

end
%% transfer task 
nstates = numel(unique(s));
% Qf = l.Q(nstates-4+1:nstates,[2,1])';
% Qf = Qf(:);
% Qf(1)=l.Q(nstates-4+1,2);
% Qf(2)=l.Q(nstates-4+1,1);
% Qf(3)=l.Q(nstates-4+2,2);
% Qf(4)=l.Q(nstates-4+2,1);
% Qf(5)=l.Q(nstates-4+3,2);
% Qf(6)=l.Q(nstates-4+3,1);
% Qf(7)=l.Q(nstates-4+4,2);
% Qf(8)=l.Q(nstates-4+4,1);

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

C = nchoosek(1:8,2);
ppc = zeros(1,8);
for k = 1: length(C)
    
    [dum_action,PC] = l.chooseActionPostLearning(Qf(C(k,2)),Qf(C(k,1))); 
        if dum_action ~=2
            PC = 1-PC; 
        end %probability of choosing 2 (here 2 =/= correct)
    ppc(C(k,1)) = ppc(C(k,1)) + PC./7;
    ppc(C(k,2)) = ppc(C(k,2)) + (1-PC)./7;
end

pll = NaN(112,1);
for i = 1:length(aa)
    dQ_post(i) = Qf(ss(i,2)) - Qf(ss(i,1));
    V_post(i) = (Vf(ss(i,2)) + Vf(ss(i,1)))/2;
    Q_post(:,i) = Qf(ss(i,:));

    if (aa(i)) && (aa(i))~=1.5                                                                          % if a choice was performed in time at the first level
        [dum_action,trl_p] = l.chooseActionPostLearning(Qf(ss(i,1)),Qf(ss(i,2))); 
        if dum_action~= aa(i)
            trl_p = 1-trl_p;
        end
        pll(i) = trl_p;          % the likelihood of the choice
    end
end

end
