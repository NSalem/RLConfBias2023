function [data_out,post]= AggregateData2(data_dir,subjects,n_sess, outcomes,reversalIdx,fileString)
if nargin<5
    reversalIdx =0;
end
if nargin<6
    fileString = 'Sub%d_Session%d.mat';
end
%%% Aggregates data over participants
%%% input
%%%     data_dir: (string) directory where to find files
%%%     subjects: (vector) which subjects to include 
%%%     n_sess: (int) number of sessions to include
%%%     outcomes: (vector) vector of length 2 containig the absolute magnitudes
%%%     of possible outcomes for these subjects. Can be [0,0.5] or [0.1,1]
%%% output
%%%     data_out: (struct) data structure containing behavioral variables
%%%     for learning task across several participants (subject x session x trial number)
%%%         Note that for choice we take missing values (lack of response)
%%%         as choice of the worst option
%%%     post: (struct) data structure containing behavioral variables of transfer task 
%%%     accross several participants (subject x trial)
    n_sub = length(subjects);
    
    for k_sub       = 1:n_sub
        
    flnm2    = strcat('PostTest',num2str(subjects(k_sub)));
    Pdata = load(fullfile(data_dir,flnm2));
    
    for k_symb = 1:8;
        S1 = Pdata.data(:,3) == k_symb;
        S2 =  Pdata.data(:,4) == k_symb;
        post.pref(k_sub,k_symb) = 100*(mean(Pdata.data(S1,6)==-1) + mean(Pdata.data(S2,6)==1))./2;
    end
    
    post.ss(k_sub,:,:)  = Pdata.data(:,3:4);
    post.aa(k_sub,:)  = Pdata.data(:,6)/2+1.5;
    post.conf(k_sub,:) = Pdata.data(:,7);
    post.rt(k_sub,:)  =Pdata.data(:,9);

        for k_sess  = 1:n_sess
            flnm    = sprintf(fileString,(subjects(k_sub)),(k_sess));
%             flnm    = strcat('Sub',num2str(subjects(k_sub)),'_Session',num2str(k_sess));
            load(fullfile(data_dir,flnm));
            
            if all(unique(outcomes) == [0,0.5])
                out_sub = double(data(:,8)/2-0.5.*(data(:,4)>2));
                cou_sub = double(data(:,9)/2-0.5.*(data(:,4)>2));
            elseif all(unique(outcomes) == [0.1,1])
                out_sub = double(data(:,4)<3).*(data(:,8) + .1.*(1-data(:,8))) - (1-double(data(:,4)<3)).*(.1*data(:,8) + (1-data(:,8)));
                cou_sub = double(data(:,4)<3).*(data(:,9) + .1.*(1-data(:,9))) - (1-double(data(:,4)<3)).*(.1*data(:,9) + (1-data(:,9)));    
            end
            if any(isnan(out_sub(:)))
            keyboard()
            end
            out_good = data(:,8);
            out_good(out_good==0)=-1;
            %%%make counterfactual nan for partial
            cou_sub(mod(data(:,4),2)>0) = NaN;
            cou_good = data(:,9);
            cou_good(mod(data(:,4),2)>0) = NaN;
            cou_good(cou_good==0)=-1;
            
            val = double(data(:,4)<3).*1 + double(data(:,4)>=3).*-1;
            info = double(mod(data(:,4),2)==0);
           
            if reversalIdx
                bestoption = zeros(size(data,1),1);
                correct = zeros(size(data,1),1);    
                for icon = [1:4]
                   idCon = find(data(:,4)==icon);
                   if mod(icon,2)==0
                       idNoRev = idCon(1:reversalIdx-1);
                       idRev = idCon(reversalIdx:end);
                       correct(idNoRev) = data(idNoRev,7);
                       correct(idRev) = 1-data(idRev,7);
                       bestoption(idNoRev) = 2;
                       bestoption(idRev) = 1;
                   else
                       bestoption(idCon) = 2;
                       correct(idCon) = data(idCon,7);
                   end
                end
            else
                correct = data(:,7);
                bestoption = ones(size(data,1),1)*2;
            end
            
            data_out.corr(k_sub,k_sess,:)    =  correct;
            data_out.conf(k_sub,k_sess,:)    =  data(:,10);
            data_out.con(k_sub,k_sess,:)     = data(:,4);
            data_out.val(k_sub,k_sess,:)     = val;
            data_out.info(k_sub,k_sess,:)     = info;
            data_out.cho(k_sub,k_sess,:)     = data(:,7)+1;
            data_out.out(k_sub,k_sess,:)     = out_sub;
            data_out.cou(k_sub,k_sess,:)     = cou_sub;
            data_out.out_good(k_sub,k_sess,:)     = out_good;
            data_out.cou_good(k_sub,k_sess,:)     = cou_good;
            data_out.RT (k_sub,k_sess,:)     = data(:,12);
            data_out.exptime(k_sub,k_sess,:)  = data(:,5);
            data_out.bestoption(k_sub,k_sess,:) = bestoption;
            data_out.duration(k_sub,k_sess) = data(end,5);
            data_out.moneycond{k_sub,k_sess} = money;
            data_out.money(k_sub,k_sess)= sum(money);
        end
    end
    data_out.dims = 'subj_sess_trial';
    [~,expname] = fileparts(data_dir);
    data_out.exp = expname;
    data_out.subjects = subjects;
end