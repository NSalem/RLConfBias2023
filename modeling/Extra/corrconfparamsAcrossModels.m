%%% correlate difficulty and Qc bias across different difficulty models for confidence %%%

regdQL = load('Results\reg_conflogit_learning.mat','confcoeff');
regdQabsL = load('Results\reg_conflogit_learning_dqabs.mat','confcoeff');
regdQGBabsL = load('Results\reg_conflogit_learning_dqgoodbad.mat','confcoeff');
regdPCL = load('Results\reg_conflogit_learning_pc.mat','confcoeff');

regdQT = load('Results\reg_conflogit_posttest.mat','confcoeff');
regdQabsT = load('Results\reg_conflogit_posttest_dqabs.mat','confcoeff');
regdQGBabsT = load('Results\reg_conflogit_posttest_dqgoodbad.mat','confcoeff');
regdPCT = load('Results\reg_conflogit_posttest_pc.mat','confcoeff');


cdiffL = corr([regdQL.confcoeff{8}(:,[2]),regdQabsL.confcoeff{8}(:,[2]),regdQGBabsL.confcoeff{8}(:,2),regdPCL.confcoeff{8}(:,[2])]);
cqcL = corr([regdQL.confcoeff{8}(:,[3]),regdQabsL.confcoeff{8}(:,[3]),regdQGBabsL.confcoeff{8}(:,3),regdPCL.confcoeff{8}(:,3)]);

cdiffT = corr([regdQT.confcoeff{2}(:,[2]),regdQabsT.confcoeff{2}(:,[2]),regdQGBabsT.confcoeff{2}(:,2),regdPCT.confcoeff{2}(:,[2])]);
cqcT = corr([regdQT.confcoeff{2}(:,[3]),regdQabsT.confcoeff{2}(:,[3]),regdQGBabsT.confcoeff{2}(:,3),regdPCT.confcoeff{2}(:,3)]);



figure()
mats  = {cdiffL,cqcL,cdiffT,cqcT};
for ip = 1:4
subplot(2,2,ip)
imagesc(mats{ip})
xticklabels({'Q_c-Q_u','|Q_c-Q_u|','Q_{GOOD}-Q_{BAD}','pc'})
xtickangle(90)
yticklabels({'Q_c-Q_u','|Q_c-Q_u|','Q_{GOOD}-Q_{BAD}','pc'})
colorbar()
end

labels = {'Q_c-Q_u','|Q_c-Q_u|','Q_{GOOD}-Q_{BAD}','pc'};

figure()
for ip = 1:4
    switch ip  
        case 1

            mymat = [regdQL.confcoeff{11}(:,[2]),regdQabsL.confcoeff{11}(:,[2]),regdQGBabsL.confcoeff{11}(:,2),regdPCL.confcoeff{11}(:,[2])];
            t = '\beta_{diff} learning';
        case 2
            mymat = [regdQL.confcoeff{11}(:,[3]),regdQabsL.confcoeff{11}(:,[3]),regdQGBabsL.confcoeff{11}(:,[3]),regdPCL.confcoeff{11}(:,[3])];
            t = '\beta_{Qc} learn';
        case 3
            t = '\beta_{diff} transfer';
            mymat =[regdQT.confcoeff{11}(:,[2]),regdQabsT.confcoeff{11}(:,[2]),regdQGBabsT.confcoeff{11}(:,2),regdPCT.confcoeff{11}(:,[2])];
        case 4
            t = '\beta_{Qc} transfer';
            mymat =[regdQT.confcoeff{11}(:,[3]),regdQabsT.confcoeff{11}(:,[3]),regdQGBabsT.confcoeff{11}(:,[3]),regdPCT.confcoeff{11}(:,[3])];
    end
    
    subplot(2,2,ip)  
    [~,ax] = plotmatrix(mymat);
    title(t)
    
    for ip1 = 1:size(ax,1)
            ylabel(ax(ip1,1),labels{ip1})
            xlabel(ax(4,ip1),labels{ip1})
        for ip2 = 1:size(ax,2)
            xlim(ax(ip1,ip2),[-10,10])
            ylim(ax(ip1,ip2),[-10,10])
        end
    end
end
