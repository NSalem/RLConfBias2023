%%%BMC family-wise comparison of confidence models
%%%compare models using previous confidence with those not using it

for itask = 1:2
    if itask ==1
        load('Results\reg_conflogit_learning_dqabs.mat','BICCONF','whichLearnModel','whichConfModel','isConfPrev')
        taskname = 'LL';
    elseif itask ==2
        load('Results\reg_conflogit_posttest_dqabs.mat','BICCONF','whichLearnModel','whichConfModel','isConfPrev')
        taskname = 'TT';
    end
    famnames ={'No conf_{t-1}','Conf_{t-1}'};
    
    [r,c] = find(~isnan(BICCONF));
    dum = BICCONF(unique(r),:);
    ind = whichLearnModel==11;
    dum = dum(:,ind);
    
    options = struct();
    options.families = {find(~isConfPrev(ind)),find(isConfPrev(ind))};
    options.DisplayWin = false;
    [postBMC,outBMC] = VBA_groupBMC(-dum'/2,options);
    BMCres{itask} = outBMC;
    figure()
    hold on
    
    bar(squeeze(100*outBMC.families.ep),...
        'FaceColor',.7*[1,1,1],...
        'EdgeColor','none')
    alpha(0.5)

    errorbar(100*outBMC.families.Ef,100*sqrt(outBMC.families.Vf),'d',...
        'Color',.5*[1,1,1],...
        'MarkerFaceColor',.7*[1,1,1],...
        'MarkerEdgeColor',.5*[1,1,1],...
        'LineStyle','None')
    plot([-1,numel(outBMC.families.Ef)+1],100*[1./numel(outBMC.families.Ef),1./numel(outBMC.families.Ef)],'r--');
    plot([-1,size(outBMC.families.Ef,1)+1],100*[0.95,0.95],'b--')
    
    set(gca,'YLim', [0 100],...
        'XTick',1:numel(outBMC.families.Ef),...
        'XTickLabels',famnames,...
        'XLim',[0 size(outBMC.families.Ef,1)+1],...
        'FontName','Arial','FontSize',14)
    xtickangle(90)
    hL = legend('Exceedence P.','Expected F.');
    hY = ylabel('Probability');
    set([hY,hL],'FontName','Arial','FontSize',14)
    set(hL,'Location','northeastoutside')
    legend boxoff  
    saveas(gcf,['Plots\BMCConfPrev',taskname,'.svg']);
    clear BICCONF whichLearnModel whichConfModel
end
