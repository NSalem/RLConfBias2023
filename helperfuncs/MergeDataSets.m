f1 = load('Results/data_all_CONF');
f2 = load('Results/data_all_sample2');

for iexp =1:numel(f2.data_all)
    f2.data_all(iexp).expN = f2.data_all(iexp).expN+numel(f1.data_all);
end


fnames = fieldnames(f1.data_all);
for ifield = 1:numel(fieldnames(f1.data_all))
    
    for iexp = 1:numel(f2.data_all)
        f2.data_all(iexp).conf = f2.data_all(iexp).corr*NaN;
        f2.data_all(iexp).conf_post = f2.data_all(iexp).aa*NaN;
    end
    
    if ~isfield(f2.data_all,fnames{ifield})
        for iexp =1:numel(f2.data_all)
            f2.data_all(iexp).(fnames{ifield}) = [];%f1.data_all(1).(ifield)*NaN
        end
    end
end


fnames = fieldnames(f2.data_all);
for ifield = 1:numel(fieldnames(f2.data_all))
    if ~isfield(f1.data_all,fnames{ifield})
        for iexp =1:numel(f1.data_all)
            f1.data_all(iexp).(fnames{ifield}) = [];%f1.data_all(1).(ifield)*NaN
        end
    end
end

f1.data_all = orderfields(f1.data_all);
f2.data_all = orderfields(f2.data_all);

data_all = [f1.data_all;f2.data_all']

save('Results/data_all.mat','data_all')