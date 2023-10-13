function CV = cvLooKsr(x,y,h,N)
    %%% cross-validation folds
    %%% x: predictor variable
    %%% y: response variable
    %%% h: smoothing window for non-parametric regression 
    %%% N: number of points for non-parametric regression 
    
    CVerr = nan(size(x));
    for i =1:numel(x)
        %%% train-test split
        xTrain = x;
        xTrain(i) = [];
        xTest = x(i);
        yTrain = y;
        yTrain(i) = [];
        yTest = y(i);
        
        %%% do kernel regression on training
        npreg = ksr(xTrain,yTrain,h,N,[min(x),max(x)]);
        
        %%% get CV error
        [~,idc] = min(abs(npreg.x-xTest));
        CVerr(i) = (yTest-npreg.f(idc))^2;
    end
    CV = mean(CVerr)';
end