function [cvMean,cvSD] = cvKfoldKsr(x,y,k,h,N)
    %%% cross-validation folds
    %%% x: predictor variable
    %%% y: response variable
    %%% k: number of folds
    %%% h: smoothing window for non-parametric regression 
    %%% N: number of points for non-parametric regression 
    
    idx = randperm(numel(x));
    xNew = x(idx);
    yNew = y(idx);
    foldSize = round(numel(x)/k);
    
    CVerr = nan([k,foldSize]);
    
    for ifold =1:k
        %%% train-test split
        startInd =(ifold-1)*foldSize+1;
        endInd =  (ifold-1)*foldSize+foldSize;
        if endInd>numel(xNew)
            endInd = numel(xNew);
        end
        idK = startInd:endInd;
        xTrain = xNew(idK);
        xTest = xNew;
        xTest(idK) = [];
        yTrain = yNew(idK);        
        yTest = yNew;
        yTest(idK) = [];
        
        %%% do kernel regression on training
        npreg = ksr(xTrain,yTrain,h,N,[min(x),max(x)]);
        
        for iobs = 1:numel(xTest)
            %%% get CV error
            [~,idc] = min(abs(npreg.x-xTest(iobs)));
            CVerr(ifold,iobs) = (yTest(iobs)-npreg.f(idc))^2;
        end
    end
    CV = mean(CVerr,2)';
    cvMean = mean(CV);
    cvSD = std(CV);
end