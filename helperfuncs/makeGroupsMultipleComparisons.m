function congroupletters = makeGroupsMultipleComparisons(c)
%%% c: output from multcompare

ncond = numel(unique(c(:,[1,2])));

letters = lower(char([1:ncond] + 64));
%       letters = {'a','b','c','d','e','f','g','h'}
%         congroupletters = {'','','','','','','',''};
congroupletters = cell(1,ncond);
for icell = 1:numel(congroupletters)
    congroupletters{icell} = '';
end

%         letters2 = ['a','b','c','d','e','f','g','h'];
prevgroups = {};
letterind = 0;
for icond = 1:ncond %iterate accross conditions
    %%% find conditions not singificantly different
    condidcs = find(any((c(:,[1,2]) == icond)'));
    groupidcs = c(condidcs,6)>.05; %%find conditions non-significantly different from ocind
    %                 groupconds = find(groupidcs);
    condsgroup = unique([c(condidcs(groupidcs),[1,2])])';
    condsgroup = condsgroup(:)';
    condsgroup = unique([condsgroup,icond]); %% this can be removed to show only groups of more than 1 condition %%%
  
    for igroupcond =condsgroup
          thisGroup = igroupcond;
        %%% for each condition non-sig different from icond,
        %%% check significance between them. If they are not,
        %%% they remain in the same group (e.g. 3 conditions
        %%% that are n.s. with each other), otherwise create
        %%% new group.
%         if igroupcond == icond
%             continue
%         end
        
        isNewSubGroup=0;
%         iminus1 = condsgroup(find(igroupcond == condsgroup)-1);
        for igroupcondrest = condsgroup%condsgroup(1:find(condsgroup ==iminus1))
            if igroupcondrest == igroupcond
                continue
            end
            %test significance with all previous. If not, use
            %new letter to indicate new group
            idcomp = find(all((c(:,[1,2]) == sort([igroupcond,igroupcondrest]))'));
            if c(idcomp,6)<=0.05
                isNewSubGroup = 1;
            else
                thisGroup = [thisGroup,igroupcondrest];
            end
        end
        if ~isNewSubGroup & ~any(thisGroup == igroupcond)
            thisGroup = [thisGroup,igroupcond];
        end
        
        if isempty(prevgroups)
            prevgroups{1} = thisGroup;
        else
            knownGroup = 0;
            for iprev = 1:numel(prevgroups)
                if (numel(prevgroups{iprev}) == numel(thisGroup) && all(sort(prevgroups{iprev}) == sort(thisGroup)))
                    knownGroup =1;
                    break
                end
            end
            if ~knownGroup 
                prevgroups{end+1} = thisGroup;
            end
            
        end
           
            
    end
        
end


    bad = zeros(numel(prevgroups),1);
    for  iprev1 = 1:numel(prevgroups)
        for iprev2 = 1:numel(prevgroups)
            if iprev2 == iprev1
                continue;
            end
            if numel(prevgroups{iprev2})~= numel(prevgroups(iprev1))
               %%% check if prev2 in prev1
               if numel(intersect(prevgroups{iprev1},prevgroups{iprev2})) == numel(prevgroups{iprev2})
                  bad(iprev1) = 1;
               end
            end
        end
    end

    newgroups = {};
    for iprev1 = 1:numel(prevgroups)
        if ~bad(iprev1)
            if isempty(newgroups)
                newgroups{1} = prevgroups{iprev1};
            else
                newgroups{end+1} = prevgroups{iprev1};
            end
        end
    end


    for igroup = 1:numel(newgroups)
        for icond = newgroups{igroup}
            congroupletters{icond} = [congroupletters{icond},letters(igroup)];
        end
    end

end