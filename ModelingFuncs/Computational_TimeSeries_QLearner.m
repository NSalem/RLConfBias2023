function [Q,V,pc,PE,ppc,pll,dQ_post,V_post,Q_post,Qf,Vf] = Computational_TimeSeries_QLearner(params,s,a,r,c,aa,ss)

%simulate hidden variables given actions and parameters, with output sorted per condition (i.e. matrices of condition x session x trial)
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


l = qLearner(params);
% nsessions=floor(length(a)/96);
% Q       = zeros(4*nsessions,2,24);        % Initial option values (all Models) as a function of conditio ("s")
% V       = zeros(4*nsessions,24);        % Initial Context values (Models 3) as a function of conditio ("s")
% pc      = zeros(4*nsessions,24);
% cc      = NaN(length(c),1);
% trialc = zeros(4*nsessions,1);                                                                % loglikelyhood
if sum(size(s)>1)==2
    ntrials = size(s,1);
else
    ntrials = 24;
end
nsessions=floor(length(s(:))/(ntrials*4));
Q       = zeros(4*nsessions,2,ntrials);        % Initial option values (all Models) as a function of conditio ("s")
V       = zeros(4*nsessions,ntrials);        % Initial Context values (Models 3) as a function of conditio ("s")
pc      = ones(4*nsessions,ntrials)*.5;
dQ_post = NaN(length(ss),1);
V_post = NaN(length(ss),1);
Q_post = NaN(2,length(ss));
trialc = zeros(4*nsessions,1);                                                                % loglikelyhood


%maybe preallocate Q spaces in learner properties...

for i = 1:length(a)
    trialc(s(i)) =  trialc(s(i)) +1;
    if ~isnan(a(i))
        %%% get probability of action being correct given current values
        [dum_action,trl_pc] = l.chooseAction(s(i));
        if dum_action ~=2
            trl_pc = 1-trl_pc;
        end
        pc(s(i),trialc(s(i))) = trl_pc;
        COMP = (1-mod(s(i),2));
        PAR = round(mod(s(i),2));
        %%% learn
        if COMP
            l.learn(s(i),a(i),r(i),c(i));
        else
            l.learn(s(i),a(i),r(i));
        end
        
        if i<length(a)
            PE(s(i),trialc(s(i))+1) = r(i)-l.V(s(i)) - l.Q(s(i),a(i)); %%%TO BE CHECKED
        end
        Q(s(i),a(i),trialc(s(i))+1) = l.Q(s(i),a(i));
        Q(s(i),3-a(i),trialc(s(i))+1) = l.Q(s(i),3-a(i));
        V(s(i),trialc(s(i))+1) = l.V(s(i));
        
    else
        Q(s(i),1,trialc(s(i))+1) = Q(s(i),1,trialc(s(i)));
        Q(s(i),2,trialc(s(i))+1) = Q(s(i),2,trialc(s(i)));
        
    end
    
end
%% now post-test
Qf(1)=Q(end-4+1,2,end); %PG75
Qf(2)=Q(end-4+1,1,end); %PG25
Qf(3)=Q(end-4+2,2,end); %CG75
Qf(4)=Q(end-4+2,1,end); %CG25
Qf(5)=Q(end-4+3,2,end); %PL25
Qf(6)=Q(end-4+3,1,end); %PL75
Qf(7)=Q(end-4+4,2,end); %CL25
Qf(8)=Q(end-4+4,1,end); %CL75


Vf(1)=V(end-4+1,end); %PG75
Vf(2)=V(end-4+1,end); %PG25
Vf(3)=V(end-4+2,end); %CG75
Vf(4)=V(end-4+2,end); %CG25
Vf(5)=V(end-4+3,end); %PL25
Vf(6)=V(end-4+3,end); %PL75
Vf(7)=V(end-4+4,end); %CL25
Vf(8)=V(end-4+4,end); %CL75


C = nchoosek(1:8,2);
ppc = zeros(1,8);
for k = 1: length(C)
    
    [dum_action,PC] = l.chooseActionPostLearning(Qf(C(k,2)),Qf(C(k,1)));
    if dum_action ~=2
        PC = 1-PC;
    end; %probability of choosing 2 (here 2 =/= correct)
    ppc(C(k,1)) = ppc(C(k,1)) + PC./7;
    ppc(C(k,2)) = ppc(C(k,2)) + (1-PC)./7;
end

pll = NaN(112,1);
for i = 1:length(aa)
    dQ_post(i) = Qf(ss(i,aa(i))) - Qf(ss(i,3-aa(i)));
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
