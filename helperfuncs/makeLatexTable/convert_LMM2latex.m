function latexCode = convert_LMM2latex(mdl,modelNames,info,round_digits,label,orientation)
%   Create a latex table code containing the output of the linear
%   mixed-effect model (created by fitlme, fitglme).

%    latexCode = convert_LMM2latexTable(mdl,modelNames,info,round_digits,orientation)

%   'latex_code'    - table latex code
%   'mdl'           - output of fitlme/fitglme. If more than one model are to
%                     be included in the table, then 'mdl should be a cell array
%                     containing each model as a separate cell.
%   'modelNames'    - cell array containing the names of the models. If left
%                     empty ([]), the default names are model 1, model 2, etc.
%   'info'          - cell array containing the list of the information
%                     required. By default, the parameter estimate, SE and
%                     p-val are included. For now the other possible arguments
%                     are 'tStat' for the t-statistic, 'F' for the
%                     F-statistic, 'n' for the number of observations, and
%                     'adjR2' for the adjusted R-squared, the name of the
%                     model criterion required for any of the model
%                     criteria.
%   'round_digits'  - number of digits included in the rounding
%   'label'         - Name used to reference the table.
%   'orientation'   - orientation of the table. By default, the orientation
%                     is vertical ('vertical'). It can be rotated by 90 deg
%                     by setting to 'rotated'.

%   Written by Pierre Petitet (pierre.petitet@gmail.com)
%   Last update: 04/07/2020

if nargin<5 label='insertLabel'; end 
if nargin<6 orientation='vertical'; end % by default, the orientation of the table is vertical
% Convert the fitglm output to a useable one for this function:
results = prepare_glme4latex(mdl,info);
nModels=length(results); % Number of models

%% Create column headers (name of the models):
colHeaders=[];
if ~isempty(modelNames) && length(unique(modelNames)) == length(modelNames) % all models have different names
    for m=1:nModels
        colHeaders=[colHeaders,'& \textbf{',modelNames{m},'}'];
    end
    
elseif ~isempty(modelNames) && length(unique(modelNames)) < length(modelNames) % some models have the same name
  uniqueModelNames = unique(modelNames);
    for m = 1:length(uniqueModelNames)
        n = sum(strcmp(modelNames,uniqueModelNames{m}));
        colHeaders=[colHeaders,'& \multicolumn{' num2str(n) '}{c}{',uniqueModelNames{m},'}'];
    end
else
    for m=1:nModels
        modelNames{m}=['Model ' num2str(m)];
        colHeaders=[colHeaders,' & \textbf{Model ',num2str(m),'}'];
    end
end
colHeaders=[colHeaders,' \\'];

%% Create row names (unique variables in the model)
all_coefNames=[];K=[];
for m=1:nModels
    all_coefNames=[all_coefNames;results{m}.coefNames];
    K=vertcat(K,repmat(m,[length(results{m}.coefNames) 1]));
end
coefNames=unique(all_coefNames);

% Re-order variables so that the ones common to all models come first:
KK = zeros([length(coefNames) nModels]);
for i = 1:length(coefNames)
    nRepet(i,1)=-sum(strcmp(all_coefNames,coefNames{i}));
    nRepet(i,2) = sum(K(strcmp(all_coefNames,coefNames{i})));
end
[~,order]=sortrows(nRepet,[1,2]); coefNames = coefNames(order);

% Names of the extra info (bottom of the table)
extraNames=[];
if isfield(results{1},'extra') extraNames=fieldnames(results{1}.extra); end

% Size of the text
if length(coefNames) <= 8 % small number of rows
    txtSize = '\footnotesize'; % large number of rows
else txtSize = '\tiny'; end

%% Create the core of the table
nLines=numel(fieldnames(results{1}.output))+1; % number of lines per cell

% First column:
for i = 1:nLines.*length(coefNames)
    if mod(i,nLines)==1
        tableCode{i}= ['\hline ', coefNames{ceil(i./nLines)}, '&'];
    else
        tableCode{i}='&';
    end
end

% Info to show that is not param estimate, SE, or p-val:
extraOutput = fieldnames(results{1}.output);
extraOutput(strcmp(extraOutput,'coef'))=[]; extraOutput(strcmp(extraOutput,'pval'))=[];

nEmptyCells=0;
% Other columns:
for m = 1:nModels
    for c = 1:size(coefNames)
        idx = find(strcmp(results{m}.coefNames,coefNames{c}));
        
        if ~isempty(idx)
            % First two lines are always beta and SE
            beta = sprintf(['%+0.' num2str(round_digits+1) 'g'],results{m}.output.coef.value(idx));
            SE = sprintf(['%0.' num2str(round_digits+1) 'g'],results{m}.output.coef.SE(idx));
            if contains(beta,'e')
                beta = [strrep(beta,'e',' \times 10^{'),'}'];
            end
            if contains(SE,'e')
                SE = [strrep(SE,'e',' \times 10^{'),'}'];
            end
            i = (c-1).*nLines+1;
            tableCode{i} = [tableCode{i},'$\beta=', beta, '$&'];
            i = (c-1).*nLines+2;
            tableCode{i} = [tableCode{i},'$SE=', SE, '$&'];
            
            
            
            % Last line is always p-value
            i = c.*nLines;
            pval = results{m}.output.pval(idx);
            if pval < 0.0001
            tableCode{i} = [tableCode{i},'\textbf{\textit{p\textless0.0001}}&'];
                
            elseif pval>=0.0001 & pval<0.05
                tableCode{i} = [tableCode{i},...
                    '\textbf{\textit{p = ',sprintf(['%0.' num2str(round_digits) 'g'],pval),'}}&'];
            else
                tableCode{i} = [tableCode{i},...
                    '$p=',sprintf(['%0.' num2str(round_digits) 'f'],pval),'$&'];
            end
            
 
            
            % Lines in-between
            for e = 1:nLines-3
                i = (c-1).*nLines+2+e;
                switch extraOutput{e}
                    case 'tStat'
                        tableCode{i} = [tableCode{i},'$',...
                            sprintf(['t_{%0.0f}=%+0.' num2str(round_digits),'f'],results{m}.output.tStat(idx,2),results{m}.output.tStat(idx,1)),'$&'];
                    case 'F'
                        tableCode{i}=[tableCode{i},'$',...
                            sprintf(['F_{%0.0f,%0.0f}=%0.',num2str(round_digits),'f'],...
                            results{m}.output.F(idx,1),results{m}.output.F(idx,2),...
                            results{m}.output.F(idx,3)),'$&'];
                end
            end
            
        else % In case this doesn't exist: empty cell
            for i = (c-1).*nLines+1 : c.*nLines
                tableCode{i}=[tableCode{i},'&'];
            end
            nEmptyCells = nEmptyCells+1;
        end
    end
end

for i = 1:nLines.*length(coefNames)
    tableCode{i}(end)=[]; % remove last '&'
    tableCode{i} = [tableCode{i},'\\']; % and replace it by a \\
end

%% Create extra lines bellow the table
tableCode_extra=[];
for i = 1:length(extraNames)
    tableCode_extra{i}=[extraNames{i},'&'];
    tableCode_extra{i}=strrep(tableCode_extra{i},'adjR2','$adj-R^2$');
    tableCode_extra{i}=strrep(tableCode_extra{i},'n','$N_{obs}$');
    
    for m=1:nModels
        if strcmp(extraNames{i},'n')
            tableCode_extra{i}=[tableCode_extra{i},...
            num2str(results{m}.extra.(extraNames{i})),'&'];
        else
        tableCode_extra{i}=[tableCode_extra{i},...
            sprintf(['%0.' num2str(round_digits) 'f'],...
            results{m}.extra.(extraNames{i})),'&'];
        end
    end
    tableCode_extra{i}(end)=[]; % remove last '&'
    tableCode_extra{i}=[tableCode_extra{i},'\\'];
end
if ~isempty(tableCode_extra) tableCode_extra{1} = ['\hline \hline ', tableCode_extra{1}]; end
%% Create title and caption
% The caption contains the formulas of the models.

formulas = {};
for m =1:numel(results)
    formulas{m} = results{m}.formula;
end

if (nEmptyCells==0 && numel(unique(formulas))==1)% means that all columns use the same model
    ttl = ['\caption{\textbf{Insert Title Here.} For each COLUMN, the model was specified as follows: ',...
        results{m}.formula, '; '];
else
    ttl = '\caption{\textbf{Insert Title Here.} Models were specified as follows. ';
    for m = 1:nModels
        ttl = [ttl, modelNames{m}, ': ', results{m}.formula '; '];
    end
end
ttl = [ttl(1:end-2) '.}'];
ttl = strrep(ttl,'|','\textbar \');
ttl = strrep(ttl,'~','$\sim$');

%% Create environment
tabvar='\begin{tabular}{r|';
for m=1:nModels
    tabvar=[tabvar,'l'];
end
tabvar=[tabvar,'}'];

%% Put together the code
if strcmp(orientation,'rotated')
    startTable = '\begin{sidewaystable}';
    endTable = '\end{sidewaystable}';
    package_required = {...
        '%   Please add the following required package to your document preamble:';...
        '%   \usepackage{rotating}'};
else
    startTable='\begin{table}';
    endTable='\end{table}';
    package_required = [];
end

latexCode=[startTable;['\centering ' txtSize];tabvar; ...
    '\hline \hline';colHeaders;'\hline';tableCode';...
    tableCode_extra' ;'\hline \hline';'\end{tabular}';ttl;['\label{tab:' label '}'];endTable];

%% Add some metadata
c=clock;
time_char=strcat(num2str(c(2)),'/',num2str(c(3)),'/',num2str(c(1)),'-',num2str(c(4)),':',num2str(c(5)),':',num2str(c(6)));
mystack = dbstack;
ThisFileNameWithPath = which(mystack(1).file);

latexCode=['%   Table generated by the matlab function: convert_LMM2latexTable.m';
    '%   Written by Pierre Petitet (pierre.petitet@gmail.com)';...
    ['%   Stored in: ' ThisFileNameWithPath]; ...
    ['%   Date created:',' ', time_char];package_required;'% ';'% ';latexCode];

end

function results=prepare_glme4latex(mdl,info)
if nargin==1 info=[];end
if iscell(mdl)
    for m = 1:length(mdl)
        results{m}=extractInfo(mdl{m},info);
    end
else
    results{1}=extractInfo(mdl,info);
end
end

function result=extractInfo(mdl,info)

%%% do satterthwaite correction for t and pvals
[b,~,stats] = fixedEffects(mdl,'DFMethod','Satterthwaite');

result.formula = char(mdl.Formula); % model formula
result.coefNames = mdl.CoefficientNames'; % coefficient names
result.output.coef.value = mdl.Coefficients.Estimate; % parameter estimates
result.output.coef.SE = mdl.Coefficients.SE; % SE
result.output.pval = stats.pValue;
% result.output.pval = mdl.Coefficients.pValue; % p-val


for i = 1:length(info)
    switch info{i}
        case 'tStat'
%             result.output.tStat = horzcat(mdl.Coefficients.tStat,... % t-stat
%                 mdl.Coefficients.DF); % DF
            result.output.tStat = horzcat(stats.tStat,... % t-stat
                stats.DF); % DF
        case 'F'
            result.output.F = horzcat(mdl.anova.DF1,... % DF1
                mdl.anova.DF2,... % DF2
                mdl.anova.FStat); % F-stat
        case 'adjR2'
            result.extra.adjR2=mdl.Rsquared.Adjusted; % ajudsted R2
        case 'n'
            result.extra.n=mdl.NumObservations; % number of observation
        case {'AIC', 'BIC', 'LogLikelihood', 'Deviance'}
            result.extra.(info{i}) = mdl.ModelCriterion.(info{i});
        otherwise
            disp([info{i} ' is an unknown parameter.'])
    end
end
end