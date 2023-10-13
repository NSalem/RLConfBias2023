%Plots to evaluate  model identifiability and parameter
%recovery for choice, using results from modelIdentChoiceFamilies.m
%(Palminteri et al 2016,
%https://doi.org/10.1016/j.tics.2017.03.011)
%Plots adapted from Correa et al 2018
%(https://doi.org/10.1523/JNEUROSCI.0457-18.2018)

% load('Results\SIMS_ResultsIdentRecov.mat')


files = {'Results\SIMS_ResultsIdentFamilies.mat'};
outnames = {'10c100c'};
famnames = {{'No SYM','SYM'},{'ABS','REL'},{'0','Q_U','Other','Last'}};

for ifile = 1
    
    load(files{ifile})
       
    for ifam = 1:numel(famnames)
        modellabels = famnames{ifam};
        
        bm = bmAll{ifam};
        pxp = epAll{ifam};
        
        nmodels = size(pxp,1);
        nsims = size(pxp,3);
        % nparams = size(Rest,1);
        mean_pxp     = squeeze(mean(pxp,3));
        n_goodcl    = bm;
        
        % Fig 1
        h1 = figure('Units', 'pixels', ...
            'Position', [400 200 1000 350]);
        set(h1,'Color',[1,1,1])
        
        for k = 1:2
            subplot(1,2,k)
            
            switch k
                case 1
                    mtp = mean_pxp;
                    lbl = 'Exceedance probability (%)';
                case 2
                    mtp = n_goodcl/nsims*100;
                    lbl = '% Best FAMILY';
            end
            
            colormap(flipud(gray))
            imagesc(flipud(mtp))
            set(gca,'XTick',1:nmodels,...
                'YTick',1:nmodels,...
                'XTickLabel',modellabels,...
                'YTickLabel',{modellabels{fliplr(1:nmodels)}},...
                'FontName','Arial');
            ylabel('Simulated model')
            xlabel('Estimated model')
            xtickangle(90)
            c = colorbar;
            c.Label.String = lbl;
        end
        
        saveas(gcf,['Plots/identFamilies_',num2str(ifam),'_',outnames{ifile},'.svg'])
    end
end