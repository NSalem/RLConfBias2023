
figure()

set(gcf,'Color',[1,1,1])
lg={'\color{green}PG_{75}','\color{green}PG_{25}','\color{cyan}CG_{75}','\color{cyan}CG_{25}',...
    '\color{red}PL_{25}','\color{red}PL_{75}','\color{magenta}CL_{25}','\color{magenta}CL_{75}'};
pair_mat = [6 5;2 1;8 7;4 3;];
newIdc = pair_mat';
lg={lg{newIdc}};
samePairMat = zeros(8,8);

optP = [0.75,0.25,0.75,0.25,0.25,0.75,0.25,0.75];
optU = [optP*1+(1-optP)*0.1]; %option utility (values in exp 1 are different but same order);
optU(5:8) = -optU(5:8);            
optU = optU(newIdc(:));

for i1 = 1:8
    for i2 = 1:8
        %%% learning pair
        if any(all(sort([i1,i2])'== [1,2;3,4;5,6;7,8]'))
            samePairMat(i1,i2) = 1;
            whichCondMat(i1,i2) = find(all(sort([i1,i2])'== [1,2;3,4;5,6;7,8]'));
        end
        
        
        %%% valence
        if any(i1== [3,4,7,8]) && any(i2==[3,4,7,8])
            valMat(i1,i2) = 1;
        elseif any(i1== [1,2,5,6]) && any(i2==[1,2,5,6])
            valMat(i1,i2) = -1;
        end
        %%%
        if any(i1== [5,6,7,8]) && any(i2==[5,6,7,8])
            infoMat(i1,i2) = 1;
        elseif any(i1== [1,2,3,4]) && any(i2==[1,2,3,4])
            infoMat(i1,i2) = -1;
        end
        
        meanU(i1,i2) = mean(optU([i1,i2]));
        diffU(i1,i2) = abs(diff(optU([i1,i2])));
        absdiffU(i1,i2) = (diff(optU([i1,i2])));

    end
end
subplot(2,3,1)
imagesc(whichCondMat)
% dum = colormap(jet);
xticks(1:8); yticks(1:8);xticklabels(lg); yticklabels(lg);colorbar;caxis([0,4]);colormap('hot');
xtickangle(90)

h = colorbar;
ylabel(h,'Learning pair',...
    'FontName','Arial')
set(gca,'FontName','Arial')
subplot(2,3,2)
imagesc(valMat);
xticks(1:8); yticks(1:8);xticklabels(lg); yticklabels(lg);colormap('hot');h = colorbar;
xtickangle(90)

ylabel(h,'Condition valence',...
    'FontName','Arial');
set(gca,'FontName','Arial')
subplot(2,3,3)
imagesc(infoMat)
xticks(1:8); yticks(1:8);xticklabels(lg); yticklabels(lg);colormap('jet');h = colorbar;
xtickangle(90)

ylabel(h,'Condition information',...
    'FontName','Arial');
set(gca,'FontName','Arial')
subplot(2,3,4)
imagesc(meanU);
xticks(1:8); yticks(1:8);xticklabels(lg); yticklabels(lg);colormap('hot');h = colorbar;
xtickangle(90)

ylabel(h,'Mean Utility',...
    'FontName','Arial')
set(gca,'FontName','Arial')
subplot(2,3,5)
imagesc(diffU);
xticks(1:8); yticks(1:8);xticklabels(lg); yticklabels(lg);colormap('hot');h = colorbar;
xtickangle(90)

ylabel(h,'U(S1) - U(S2)',...
    'FontName','Arial') 
set(gca,'FontName','Arial')
subplot(2,3,6)
imagesc(absdiffU);
xticks(1:8); yticks(1:8);xticklabels(lg); yticklabels(lg);colormap('hot');h = colorbar;
xtickangle(90)
ylabel(h,'|U(S1) - U(S2)|',...
    'FontName','Arial')
set(gca,'FontName','Arial')

set(gcf,'Position',[0,0,1200,500])