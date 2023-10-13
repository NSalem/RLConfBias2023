    function modelFitPARALLEL(datafilename, outfilenameprefix, modelsinfoscript, whichexp, whichmodel, whichsess, useposttest,useParallel)
    %%% Fit all (selected) models to data from all experiments
    %%% Load data file, perform model fitting, save parameters and model evidence 
    %%% Inputs:
%         datafilename: filename of .mat data to be fitted
%         outfilenameprefix: prefix for output file 
%         modelsinfoscript: filename of script creating array of structs
%         with specifications for each model
%         whichexp: vector indicating which experiments from data file to
%         fit (default all)
%         wichmodel: vector indicating which models from modelsinfoscript
%         to fit (default all)
%         whichsess: vector indicating which experimental sessions to use
%         (default all)
%         usepost: bool, whether to use the post-test (transfer task data)
%         useparallel: bool, whether to use parallel computing for faster
%         fitting

    addpath('ModelingFuncs\');
    resultsdir = ['Results', filesep];

    if nargin<1 || isempty (datafilename)
        [file,path] = uigetfile('*.mat','Select data file');
        datafilename = [path,file];
    end
    if nargin<2 || isempty(outfilenameprefix)
        outfilenameprefix = '';
    end

    exps = load(datafilename);
    exps = exps.data_all;
    
    if nargin<3 || isempty (modelsinfoscript)
        [file,path] = uigetfile('*.m','Select model info script');
        modelsinfoscript = [path,file];
    end
     
    run(modelsinfoscript); %% load model specifications %%

    if nargin<4 || isempty(whichexp)
        whichexp = 1:numel(exps);
    end
        
    if nargin<5 || isempty(whichmodel)
        whichmodel = 1:numel(modelsinfo);
    end
    if nargin<6
        whichsess = [];
    end
    if nargin<7
        useposttest = 1;
    end
    
    if useParallel
        if isempty(gcp('nocreate'))
            parpool
        end
        parArg = Inf;
    else
        parArg = 0;
    end

    for iexp = whichexp
        disp(['Fitting exp', num2str(iexp)])
        if isempty(whichsess)
            n_sess = size(exps(iexp).con,2);
            sessthisexp = 1:n_sess;
        else
            sessthisexp = whichsess;
        end
        n_trials = size(exps(iexp).con,3);
        n_sub = size(exps(iexp).con,1);
        
            con  =  cell(n_sub,1);   % "con": conditions: 1: reward partial, 2 reward complete, 3, punishment partial, 4, punishment complete
            out  = cell(n_sub,1);
            cho = cell(n_sub,1);     % "cho": choice; 2 for correct option (most rewarding or less punishing), 1 for incorrect (less rewarding, most punishing)
            out = cell(n_sub,1);     % "out": outcome (50c, 0c or -50c of euros)
            cou = cell(n_sub,1);     % "cou": counterfactual outcome (50c, 0c or -50c of euros, not shown when con = 1 or 3)

%             ss_sub = mat2cell(NaN(n_sub,1,));
            ss  = NaN(size(exps(iexp).ss));
            aa  = NaN(size(exps(iexp).aa));
            
        for k_sub       = 1:n_sub
            for i_sess  = 1:numel(sessthisexp)
                k_sess = sessthisexp(i_sess);
                k_in = (i_sess-1)*n_trials+1;
                k_out =  i_sess*n_trials;
                con_sub(k_in:k_out,1) = (k_sess-1)*4 + exps(iexp).con(k_sub,k_sess,:);
                cho_sub(k_in:k_out,1) = exps(iexp).cho(k_sub,k_sess,:);
                out_sub(k_in:k_out,1) = exps(iexp).out(k_sub,k_sess,:);
                cou_sub(k_in:k_out,1) = exps(iexp).cou(k_sub,k_sess,:);
            end

            con{k_sub}  =  con_sub;     % "con": conditions: 1: reward partial, 2 reward complete, 3, punishment partial, 4, punishment complete
            cho{k_sub}  =  cho_sub;     % "cho": choice; 2 for correct option (most rewarding or less punishing), 1 for incorrect (less rewarding, most punishing)
            out{k_sub}  =  out_sub;     % "out": outcome (50c, 0c or -50c of euros)
            cou{k_sub}  =  cou_sub;     % "cou": counterfactual outcome (50c, 0c or -50c of euros, not shown when con = 1 or 3)

            ss(k_sub,:,1:2)  = exps(iexp).ss(k_sub,:,:);
            aa(k_sub,:)      = exps(iexp).aa(k_sub,:,:);
            
            clear con_sub cho_sub out_sub cou_sub
        end

        %% Parameters optimization
        % This part requires the Matlab Optimization toolbox
        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIter', 10000); % These increase the number of iterations to ensure the convergence
        for model = whichmodel;
            lb = modelsinfo{model}.lb;
            ub = modelsinfo{model}.ub;
            x0 = modelsinfo{model}.x0;

            dummodelsinfo = modelsinfo;
            parfor k_sub = 1:n_sub    
%             for k_sub = 1:n_sub
                if useposttest
                    aa_this_sub = aa(k_sub,:);
                    ss_this_sub = squeeze(ss(k_sub,:,1:2));
                else
                    aa_this_sub = [];
                    ss_this_sub = [];
                end
                
                % LPP (Laplace appriximation of the posterior probability)
                % optimization, minimizing the negative Laplace-Approximated Model
                % Evidence
               thismodelinfo = dummodelsinfo{model};
               [parametersLPP{k_sub,model},LPP(k_sub,model),~,reportLPP(k_sub,model),~,gradientLPP{k_sub,model},hessianLPP{k_sub,model}]=fmincon(@(x) GetModelLL_QLearner(x,thismodelinfo,con{k_sub},cho{k_sub},out{k_sub},cou{k_sub},aa_this_sub,ss_this_sub,1),x0,[],[],[],[],lb,ub,[],options);
                    thisH = hessianLPP{k_sub,model};
                    thisLPP = LPP(k_sub,model);
                    k = numel(dummodelsinfo{model}.paramnames);
                    %NEGATIVE model evidence
                    LAME(k_sub,model) =  thisLPP - k/2*log(2*pi) + real(log(det(thisH))/2);%Laplace-approximated model evidence
                    ll(k_sub,model) = GetModelLL_QLearner(parametersLPP{k_sub,model},dummodelsinfo{model},con{k_sub},cho{k_sub},out{k_sub},cou{k_sub},aa_this_sub,ss_this_sub,0)
            end
        end
    
            datainfo.experiment = exps(iexp).exp;
            datainfo.subjects{iexp} = exps(iexp).subjects;
            datainfo.sessions = sessthisexp;
            datainfo.postlearning = useposttest;
         
        fitinfo.date = date;
        save(['Results',filesep,'RLMODEL',outfilenameprefix, 'exp',num2str(exps(iexp).expN)],'parametersLPP','LAME','modelsinfo','datainfo','ll', 'fitinfo');        
        clear LAME LPP hessianLPP ll con cho out cou aa ss parameters parametersLPP modelinfo fmin_info clear ub lb x0 n_trials n_sub n_sess reportLPP gradientLPP 
    end
end