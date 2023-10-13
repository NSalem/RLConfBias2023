resultsdir = ['Results',filesep];
outfilename = 'data_all_sample2';

whichexp = 1:6;
exps(1).exp    = 'Stu1';
exps(1).subjects    = [1:28];           % 3 sess
exps(1).n_sess      = 4;
exps(1).expN        = 1;
exps(1).filenameprefix = 'TestSub';

exps(2).exp    = 'Stu2';
exps(2).subjects    = [1:15,41:49];           % 3 sess
exps(2).n_sess      = 2;
exps(2).expN        = 2;
exps(2).filenameprefix = 'Learning';

exps(3).exp    = 'Stu3';
exps(3).subjects    = [17:19,21,23,25:29,31:33,64:68,71,72];           % 3 sess
exps(3).n_sess      = 3;
exps(3).expN        = 3;
exps(3).filenameprefix = 'Sub';

exps(4).exp    = 'Stu4a';
exps(4).subjects    = [1:2:39];           % 3 sess
exps(4).n_sess      = 2;
exps(4).expN        = 4;
exps(4).filenameprefix = 'Sub';

exps(5).exp    = 'Stu4b';
exps(5).subjects    = [30:54];           % 3 sess
exps(5).n_sess      = 3;
exps(5).expN        = 5;
exps(5).filenameprefix = 'Sub';

% exps(6).exp    = 'Experiment 1';
% exps(6).subjects    = [17:20 22:35];           % 3 sess
% exps(6).n_sess      = 3;
% exps(6).expN        = 6;
% exps(6).filenameprefix = 'Sub';



for iexp =1:numel(exps)
    
    data_dir    = ['Data',filesep,exps(iexp).exp];
    subjects    = exps(iexp).subjects;           % 3 sess
	n_sub       = length(subjects);
    n_sess      = exps(iexp).n_sess;
    expN        = exps(iexp).expN;

   
    for k_sub       = 1:n_sub
        for k_sess  = 1:n_sess
            flnm    = strcat(exps(iexp).filenameprefix,num2str(subjects(k_sub)),'_Session',num2str(k_sess));
            load(fullfile(data_dir,flnm));

            k_in = (k_sess-1)*size(data)+1;
            k_out =  k_sess*size(data);

            con_sub(k_in:k_out,1) = (k_sess-1)*4 + data(:,4);
            if strcmp(exps(iexp).exp,'Stu1')
                colcho = 8;colout = 9;colcou = 10; rtcol = 11;
            else
                colcho = 7;colout = 8;colcou = 9; rtcol = 10;
            end
            exps(iexp).con(k_sub,k_sess,:)     = data(:,4);
            
            exps(iexp).val(k_sub,k_sess,:)  =  double(data(:,4)<3).*1 + double(data(:,4)>=3).*-1;
            exps(iexp).info(k_sub,k_sess,:)  =  double(mod(data(:,4),2)==0);
            
            out_sub = double(data(:,colout)/2-0.5.*(data(:,4)>2));
            cou_sub = double(data(:,colcou)/2-0.5.*(data(:,4)>2));
            
            cou_sub(mod(data(:,4),2)>0) = NaN;
                        
            exps(iexp).corr(k_sub,k_sess,:)    =  data(:,colcho);
            exps(iexp).cho(k_sub,k_sess,:)     = data(:,colcho)+1;
            exps(iexp).out(k_sub,k_sess,:)     = out_sub;
            exps(iexp).cou(k_sub,k_sess,:)     = cou_sub;
            exps(iexp).RT(k_sub,k_sess,:)     = data(:,rtcol);
            exps(iexp).out_good(k_sub,k_sess,:)     = data(:,colout);
            exps(iexp).cou_good(k_sub,k_sess,:)     = data(:,colcou);
            exps(iexp).bestoption(k_sub,k_sess,:) = ones(size(data,1),1)*2;
        end

        try
            flnm    = ['PostTest',num2str(subjects(k_sub))];
            load(fullfile(data_dir,flnm)); 
        catch
            try 
                flnm    = ['PostTest',num2str(subjects(k_sub)),'_1'];
                load(fullfile(data_dir,flnm));
            catch
                flnm    = ['PostTraining',num2str(subjects(k_sub)),'_1'];
                load(fullfile(data_dir,flnm)); 
            end
        end
        exps(iexp).ss(k_sub,:,1:2)  = data(:,3:4);
        exps(iexp).aa(k_sub,:)      = data(:,6)/2+1.5;

    end
end

data_all = exps;
save([resultsdir,outfilename],'data_all');


    
