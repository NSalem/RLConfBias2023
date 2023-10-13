% Aggregate data of 5 experiments (each 18 participants x 3 sessions)
% Resuls in a cell array of 5 elements, each containing a data structure
% with behavioral variables for one experiment

datadir = ['Data',filesep];
resultsdir = ['Results', filesep];
outfilename = 'data_all_CONF';
expdir    = 'Experiment 1';
subjects    = [17:20 22:35];           % 3 sess
n_sess      = 3;
[data1, post1] = AggregateData2([datadir,expdir],subjects,n_sess,[0,0.5]);
data1.conf = .5 + data1.conf./20;
post1.conf = .5 + post1.conf./20;
data1.expN = 1;
%exp2
expdir    = 'Experiment 2';
subjects    = [1:18];           % 3 sess
n_sess      = 3;
[data2, post2] =AggregateData2([datadir,expdir],subjects,n_sess,[0.1,1]);
data2.expN = 2;

% Experiment 6
expdir    = 'Experiment 6 CCT1';
subjects    = [1:18];           % 3 sess
n_sess      = 3;
[data6, post6] =AggregateData2([datadir,expdir],subjects,n_sess,[0.1,1]);
data6.expN = 3;

%Experiment 7
expdir    = 'Experiment 7 CCT2';
subjects    = [201:215 217:219];           % 3 sess
n_sess      = 3;
[data7, post7] =AggregateData2([datadir,expdir],subjects,n_sess,[0.1,1]);
data7.expN = 4;

%Experiment 8
expdir    = 'Experiment 8 CCT3';
%subjects    = [301:304 306:313 316:318];           
subjects    = [301:318];% 3 sess
n_sub       = length(subjects);
n_sess      = 3;
[data8, post8] =AggregateData2([datadir,expdir],subjects,n_sess,[0.1,1]);
data8.expN = 5;

data_all = [data1;data2;data6;data7;data8];
post_all = [post1;post2;post6;post7;post8];

    for iexp = 1:numel(data_all)
        data_all(iexp).conf_post = post_all(iexp).conf;
        data_all(iexp).rt_post = post_all(iexp).rt;
        data_all(iexp).aa  = post_all(iexp).aa;
        data_all(iexp).ss = post_all(iexp).ss;
    end

    save([resultsdir,outfilename],'data_all');


