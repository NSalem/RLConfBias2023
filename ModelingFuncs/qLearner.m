classdef    qLearner < handle
    properties
        %%% parameters and specifications %%%
        startQ %starting Q values
        beta %inverse temperature for choice
        betaTT %inverse temperature for choice in transfer task (default to same as beta) 
        lr1 %factual or confirmatory leraning rate (depending of specifications)
        lr2 %counterfactual or disconfirmatory learning rate
        lr3 %contextual learning rate
        lossScale = 1; %parameter for loss aversion (unused by default) 
        lossScaleType % type of loss aversion (unused by default)
        utilityExp = 1; %another parameter for loss aversion 
        wcouType = 1; %1: Rv = (Rc*(1-w)+Ru*w)/2, 2: Rv = Rc+Ru/(1+w)
        contextual %whether 
        confirmatory
        counterfactual
        RuHatType = 'None'; 
        fakeRu1
        fakeRu2
        ROtherWeight = 1;
        useRuCF = 0;
        fi = 0; % %choice trace decision bias (-inf to +inf) %see Palminteri (2023) https://doi.org/10.1037/bne0000541 
        tau =0 ; % choice trace accumulation, aka impulsivity (0 to 1)
        %         omega = 0; %normalization (0 to 1)
        
        %%% hidden variables %%%
        V %current context values
        Q %current action values
        C  %current choice trace
        s = []; %past states
        r = []; %past rewards
        c = []; %past choices
        a = []; %past actions   
    end
    
    methods
        function obj = qLearner(params)
            %Object = qLearner(params)
            %Create a q-learner object.
            %params is a struct with following fields:
%                 beta: inverse choice temperature
%                 betaTT (optional): inverse choice temperature for transfer task (by
%                 default same as beta)
%                 lr1: learning rate for confirmatory/factual outcomes
%                 lr2: (optional) learning rate for disconfirmatory/counfertactual outcomes
%                 lr3: (optional) learning rate for context value
%                 context.
%                 lossScale (optional): parameter for scaling of loss outcomes 
%                 respective to gains (unused by default) 
%                 lossScaleType (optional): type of loss scaling
%                 ('lambda','PVL','EVL', or empty, as in Lejarrega and
%                 Hertwig 2016).Unused by default
%                 utilityExp (optional): exponent of utility function (default 1)
%                 in contextual learning (default 1).
%                 contextual: bool, whether learner should learn context value and
%                 substract it from RPE (default 0)
%                 confirmatory: bool, whether to use different learning rates for
%                 confirmatory and disconformatory outcomes (positve and negative
%                 prediction errors respectively for factual option, inverted for
%                 counterfactual) (default 0).
%                 counterfactual: bool, whether using different learning rates for
%                 factual and counterfactual options. Only applied if confirmatory =
%                 false (default 0)
%                 RuHatType type of imagined counterfactual ('None','Qu', 
%                 'fakeRu', 'RMin', 'ROther', 'RLast','ROther2')
%                       'None': No imagined counterfactual
%                       'Qu': learnt value of unchosen option (as in Palminteri 2015) 
%                       'fakeRu': arbitrarily specified fictive imagined
%                       outcomes
%                       'RMin': worst experienced outcome in state
%                       'ROther': Other possible outcome in state
%                       'ROther2': Other possible outcome in state, but
%                       with "normalization" (divide by max magnitude in
%                       exp)
%                       'RLast': Last observed outcome for unchosen option
%                 fakeRu1: arbitrary value for imagined forgone outcome for
%                 reward conditions. Only used if RuHatType='fakeRu'  
%                 fakeRu2:arbitrary value for imagined forgone outcome for
%                 punishment conditions.Only used if RuHatType='fakeRu'
%                 ROtherWeight: coefficient scaling RuHat (w)
%                 fi = 0;   %choice trace decision bias (-inf to +inf) %see Palminteri (2023) https://doi.org/10.1037/bne0000541 
%                 tau =0 ; % choice trace accumulation, aka impulsivity (0 to 1)
            
            %Methods:
            %   learn
            %   chooseAction
            %   getActionLikelihood
            
            dum = properties(obj);
            for fn = fieldnames(params)'
                if any(strcmp(dum,fn))
                    obj.(fn{1}) = params.(fn{1});
                end
            end
            
            if ~isfield(params,'startQ') || isempty(params.startQ)
                obj.startQ = [0,0];
            end
            
            if ~isfield(params,'contextual')
                obj.contextual = false;
            end
            
            if ~isfield(params,'confirmatory')
                obj.confirmatory = false;
            end
            
            if ~isfield(params,'counterfactual')
                obj.counterfactual = false;
            end
            
            if ~isfield(params,'betaTT')|| isempty(params.betaTT)
                obj.betaTT = obj.beta;
            end
        end
        
        function [Qc,V] = learn(obj, s, a, r,c)
            %%% updates value of the chosen (and
            %%% unchosen option if c is given), as well as the context
            %%% value if obj.relative is true
            %%% s: integer indicating current state
            %%% a: integer (1 or 2) indicating the action taking this
            %%% trial
            %%% r: (float) factual monetary reward
            %%% c: (float, empty or NaN) counterfactual monetary reward (optional)

            %define whether feedback partial or complete (showing counterfactual outcome)
            if nargin<5 
                c = 0;
                COMP = 0;
                PAR = 1;
            else
                COMP = 1;
                PAR = 0;
            end
            
            if s>size(obj.Q,1)%||a>size(obj.Q,2)
                obj.Q(s, :) = obj.startQ; %start Q-values for new state
                obj.C(s,:) = [0,0]; %start choice trace for new state
            end
            if s>size(obj.V,2)
                obj.V(s) = 0; %start V-values
            end
            
            if obj.fi~=0 % update choice trace if using perseveration 
                obj.updateChoiceTrace(s,a) %see Palminteri (2023) https://doi.org/10.1037/bne0000541 
            end
            %% Generate imagined forgone outcome 
            if strcmp(obj.RuHatType,'fakeRu')
                if isempty(obj.fakeRu2)
                    obj.fakeRu2 = obj.fakeRu1;
                end
                ruHat = obj.fakeRu1*(mod(s,4)+4*(mod(s,4)==0)<3)+obj.fakeRu2*(mod(s,4)+4*(mod(s,4)==0)>2);
            elseif strcmp(obj.RuHatType,'RMin')
                if ~isempty(obj.s) && any(obj.s == s)
                    pastr_state = obj.r(obj.s == s);
                    pastc_state = obj.c(obj.s ==s);
                    ruHat = nanmin([pastr_state,pastc_state]);
                else
                    ruHat = 0;
                end
            elseif strcmp(obj.RuHatType,'ROther')
                if obj.V(s) ~= 0
                    if r == 1
                        ruHat = .1;
                    elseif r ==.1
                        ruHat = 1;
                    elseif r == -1
                        ruHat = -.1;
                    elseif r == -.1
                        ruHat = -1;
                    elseif abs(r) == 0.5
                        ruHat = 0;
                    elseif r ==0
                        ruHat = .5*sign(obj.V(s));
                    end
                    ruHat = ruHat.*obj.ROtherWeight; 
                else
                    ruHat = r;
                end
                
            elseif strcmp(obj.RuHatType,'ROther2') %"pseudo-normalized" outcomes
                if obj.V(s) ~= 0
                    if r == 1
                        ruHat = .1;
                    elseif r ==.1
                        ruHat = 1;
                    elseif r == -1
                        ruHat = -.1;
                    elseif r == -.1
                        ruHat = -1;
                    elseif abs(r) == 0.5
                        ruHat = 0;
                    elseif r ==0
                        ruHat = .5*sign(obj.V(s));
                    end
                else
                    ruHat = r;
                end
                
                %convert 0.5 to 1 for exp1
                if abs(r) == 0.5 || r ==0
                    r = sign(r);
                    c = sign(c);
                    ruHat = sign(ruHat);
                end
               
                ruHat = ruHat.*obj.ROtherWeight;
                
            elseif strcmp(obj.RuHatType, 'Qu')
%                 ruHat = scaleLossAverse(obj,obj.Q(s,3-a));
                ruHat = obj.Q(s,3-a);
            elseif strcmp(obj.RuHatType, 'RLast')
                whenThisChosen = find(obj.s == s & obj.a == a);
                if ~isempty(whenThisChosen)
                    lastC = whenThisChosen(end);
                    ruHat = obj.r(lastC);
                else
                    ruHat = 0;
                end
                
                %%% expand list of past states, rewards and actions
                obj.s = [obj.s s];
                obj.r = [obj.r r];
                obj.a = [obj.a a];
                if nargin>4
                    obj.c = [obj.c c];
                else
                    obj.c = [obj.c NaN];
                end
            
            elseif strcmp(obj.RuHatType,'V')
                ruHat = obj.V(s);
            elseif strcmp(obj.RuHatType,'V+Qu')
                ruHat = obj.V(s)+obj.Q(s,3-a);

            elseif strcmp(obj.RuHatType,'None')
                ruHat = 0;
            end

            %% Apply value scaling (i.e. loss aversion). Unused
            if ~isempty(obj.lossScaleType)
                ur = scaleLossAverse(obj,r);
                uc = scaleLossAverse(obj,c);
                ruHat = scaleLossAverse(obj,ruHat);
            else
                ur = r;
                uc = c;
            end
            
            %% Apply normalization (UNUSED)
%             if obj.omega >0 && obj.V(s)~=0 
%                 urREL = ur./abs(obj.V(s))+max(0,-obj.V(s)/abs(obj.V(s)));
%                 ucREL  = uc./abs(obj.V(s))+max(0,-obj.V(s)/abs(obj.V(s)));
%                 ruHatREL = ruHat./abs(obj.V(s))+max(0,-obj.V(s)/abs(obj.V(s)));
%              
%                 ur = obj.omega*urREL+(1-obj.omega)*ur;
%                 uc = obj.omega*ucREL+(1-obj.omega)*uc;
%                 ruHat = obj.omega*ruHatREL+(1-obj.omega)*ruHat;
%             end

            
            %% Factual and counterfactual RPE
            if PAR
                uc = ruHat;
            end
            deltaI =  ur - obj.V(s) - obj.Q(s,a);                                        % the prediction error for the factual choice
            deltaC =  (uc - obj.V(s) - obj.Q(s,3-a)) * (COMP);                    % the prediction error for counterfactual choice[(1-mod(s(i),2))=0 if s=1,3]
            
            %% Udapte chosen and unchosen options depending on specifications
            if obj.confirmatory %confirmatory/disconfirmatory update of factual and counterfactual options
                obj.Q(s,a) = obj.Q(s,a) + obj.lr1.*deltaI.*double(deltaI>0) + obj.lr2.*deltaI.*double(deltaI<0);                                     % the delta rule for the factual choice
                obj.Q(s,3-a) = obj.Q(s,3-a) + obj.lr2.*deltaC.*double(deltaC>0)+ obj.lr1.*deltaC.*double(deltaC<0);                                 % the delta rule for the counterfactual choice
                
            elseif obj.counterfactual %factual/counterfactual update
                obj.Q(s,a) = obj.Q(s,a) + obj.lr1.*deltaI;                                     % the delta rule for the factual choice
                obj.Q(s,3-a) = obj.Q(s,3-a) + obj.lr2.*deltaC;                                 % the delta rule for the counterfactual choice
            else % just factual
                obj.Q(s,a) = obj.Q(s,a) + obj.lr1.*deltaI;                                     % the delta rule for the factual choice
            end
            
            %% update V
            if obj.contextual && obj.wcouType ==1
                if PAR
                    deltaV =  (ur+ruHat)/2 -obj.V(s);
                else
                    deltaV =  (ur+uc)/2 -obj.V(s);
                end
                obj.V(s) = obj.V(s) + obj.lr3 .* deltaV;
            
            elseif obj.contextual && obj.wcouType ==2
                if PAR
                    Rv = (ur+ruHat*obj.ROtherWeight)/(1+obj.ROtherWeight);
                    deltaV =  Rv -obj.V(s);
                else
                    deltaV =  (ur+uc)/2 -obj.V(s);
                end
                obj.V(s) = obj.V(s) + obj.lr3 .* deltaV;
            end
            
            Qc = obj.Q(s,a); V = obj.V(s);
        end        
        
        function updateChoiceTrace(obj,s,a)
            CPE_c = 1-obj.C(s,a);
            CPE_u = 0-obj.C(s,3-a);
            obj.C(s,a) = obj.C(s,a) + obj.tau*CPE_c;
            obj.C(s,3-a) = obj.C(s,3-a) + obj.tau*CPE_u;
        end
            
        function [action,p] = chooseAction(obj, state)
            %state: number representing state
            %Outputs:
            %action: (1 or 2)
            %p: Probability of having chosen action based on the state      
            if state>size(obj.Q,1)
                obj.Q(state,:) = obj.startQ;
                obj.C(state,:) = [0,0];
            end  
            
            if state>size(obj.C,1)
                obj.C(state,:) = [0,0];
            end
            
            Q2 = obj.Q(state,2);
            Q1 = obj.Q(state,1);
            dQ = Q2 - Q1; % correct vs incorrect
            dC = obj.C(state,2)-obj.C(state,1);
%             pc = 1./(1+exp(-dQ.*obj.beta));
            pc = 1./(1+exp(-obj.beta*dQ-obj.fi*dC));
            action = double(rand<pc) + 1;
            p = pc.*(action==2) + (1-pc).*(action==1);
        end
        
        function lik = getActionLikelihood(obj,s,a,Q)
            %%%% get likeligood of given action for current values
            %%%%input: 
            %%%%s (int): state 
            %%%%a (int): action (1 or 2)
            %%%%Q (vector length 2): [Q (unchosen), Q (chosen)] 
            %%%%output:
            %%%%p (float): probability
            
           if s>size(obj.Q,1)
                obj.Q(s,:) = obj.startQ;
                obj.C(s,:) = [0,0];
                obj.V(s) = 0;
           end  
           
            if nargin<4
                Qc = obj.Q(s,a);
                Qu = obj.Q(s,3-a);
                Cc = obj.C(s,a);
                Cu = obj.C(s,3-a);
                dC = obj.C(s,a)-obj.C(s,3-a);
            else 
                Qc = Q(2);
                Qu = Q(1);
                Cc = obj.C(2);
                Cu = obj.C(1);
                dC = obj.C(s,2)-obj.C(s,1);
            end
            dQ = Qc-Qu;
            
            p = 1./(1+exp(-obj.beta*dQ-obj.fi*dC));    
            lik = log(p); 
            
        end
        
        
    function lik = getActionLikelihoodTransfer(obj,dQ)
            %%%% dQ = Q chosen - Q unchosen 
            p = 1./(1+exp(-obj.betaTT*dQ));    
            lik = log(p); 
    end
        
                
        function [action,p] = chooseActionPostLearning(obj,Q1,Q2)
            %output an action (and the current probability of having
            %chosen) based on given options
            
            if nargin<2 || isempty(Q2)
                Q2 = 0;
            end
            dQ = Q2-Q1; 
            
            pc = 1./(1+exp(-obj.betaTT.*dQ));
            action = double(rand<pc) + 1;
            
            if action == 2
                p = pc;
            else
                p = 1-pc;
            end
        end
        
        function u = scaleLossAverse(obj,x)
            if strcmp(obj.lossScaleType,'lambda')
                u = x*(x>0)+obj.lossScale*x*(x<0);
            elseif strcmp(obj.lossScaleType,'PVL')
                u = (x*(x>0))^obj.utilityExp + obj.lossScale*-abs(x*(x<0))^obj.utilityExp;
            elseif strcmp(obj.lossScaleType,'EVL')
                u = (2-obj.lossScale)*x*(x>=0)+obj.lossScale*x*(x<0);
            else
                u = x;
            end
        end
    end
    methods (Static)
    end
end
