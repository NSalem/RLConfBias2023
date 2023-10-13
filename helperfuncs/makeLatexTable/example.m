% This script demonstrates how to use 'convert_LMM2latex.m' and 'convert_regression2latex.m'
% Written by Pierre Petitet (pierre.petitet@gmail.com)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fitlm/fitglm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a dummy design matrix
designMatrix = table;
for i = 1:6
designMatrix.(['var' num2str(i)]) = rand(2000,1);
end

%% Create a dummy analysis
lm{1} = fitlm(designMatrix,'var1 ~ 1 + var2*var3');
lm{2} = fitlm(designMatrix,'var1 ~ 1 + var2*var3 + var4');
lm{3} = fitlm(designMatrix,'var1 ~ 1 + var3 + var4 + var5 + var6');

%% Create a LATEX table for a mixed-effect model
latexCode = convert_LMM2latex(...
    lm,... % array of glme results
    {},... % name of the models
    {'n','adjR2','BIC'},... % information required
    2,... % round digits
    'test_lm',... % table label
    'vertical'); % orientation (vertical or rotated)

% In case you want to replace the name of certain variables:
latexCode=replaceWords(latexCode,...
    {'Insert Title Here.','var1','var2','var3'},... % old names
    {'My amazing results.','error','trial','difficulty'}); % % new names

% Save the table
fid=fopen('test_lm.txt','w'); fprintf(fid,'%s\n',latexCode{:}); fclose(fid);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fitlme/fitglme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a dummy design matrix
designMatrix = table;
designMatrix.response = rand(2000,1);
for i = 1:6
designMatrix.(['expVariable' num2str(i)]) = rand(2000,1);
end
designMatrix.subject = repmat([1:20]',[100,1]);

%% Create a dummy analysis
for ex = 1:4 % experiment loop
lmm{ex} = fitlme(designMatrix,['response ~ 1 + expVariable5 + expVariable6 + expVariable' num2str(ex),...
    '+ (1 + expVariable5 + expVariable6 + expVariable' num2str(ex) '| subject)']);
end

%% Create a LATEX table for a mixed-effect model
latexCode = convert_LMM2latex(...
    lmm,... % array of glme results
    {'Exp. 1', 'Exp. 2', 'Exp. 3', 'Exp. 4'},... % name of the models
    {'F','n','adjR2','BIC'},... % information required
    2,... % round digits
    'test_lmm',... % table label
    'vertical'); % orientation (vertical or rotated)

% In case you want to replace the name of certain variables:
latexCode=replaceWords(latexCode,...
    {'Insert Title Here.','response','expVariable1','expVariable2'},... % old names
    {'My amazing results.','error','trial','difficulty'}); % % new names

% Save the table
fid=fopen('test_lmm.txt','w'); fprintf(fid,'%s\n',latexCode{:}); fclose(fid);