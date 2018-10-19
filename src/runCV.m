function results  = runCV(Xdat, Ydat, f_model, x0, lb, ub, opts, ilb, iub, psize, splits)
%runCV Compute performance on test set for all CV splits
%   calls ml_ga_fmincon()
%   uses objective function f_obj_adapt()

numCrossVal = size(splits,2);
numPredictors = length(x0);

% Cross-validation
sampleWeightAll = zeros(numCrossVal,numPredictors);
Rsq = zeros(numCrossVal,2);
LLH = zeros(numCrossVal,2);
RMSE = zeros(numCrossVal,2);
exitflag = zeros(numCrossVal,1);

for s=1:numCrossVal
    
    IDX = splits(:,s);
    
    % Regressors and target selection
    Xtrain = Xdat(IDX,:);
    ytrain = Ydat(IDX);
    
    % Test set
    testX = Xdat(~IDX,:);
    testY = Ydat(~IDX);
    
    % Optimization
    [ fval, ~, w, exitflag(s)] = ml_ga_fmincon( Xtrain, ytrain, f_model, x0, lb, ub, opts, ilb, iub, psize );
    
    % Model output for fitted parameters
    YmodelTrain = f_model(w,Xtrain);
    YmodelTest = f_model(w,testX);
    
    % Residuals
    trainResiduals = YmodelTrain - ytrain;
    testResiduals = YmodelTest - testY;
    
    % Estimate of the SD of the response distribution on
    % training set
    [ ll, S ] = f_obj_adapt( YmodelTrain, ytrain );
    
    % Compute log likelihood on training and test set
    LLH(s,:) = [ f_obj_adapt(YmodelTrain,ytrain,S) ...      % Log-likelihood
        f_obj_adapt(YmodelTest,testY,S) ];
    % Other measures
    RMSE(s,:) = [sqrt(mean(trainResiduals.^2)) sqrt(mean(testResiduals.^2))];       % RMSE comparison
    Rsq(s,:) = [rsquared( YmodelTrain, ytrain ) rsquared( YmodelTest, testY )];      % R^2
    
    sampleWeightAll(s,:) = w;
end

% Save data of optimization procedure
results.model = f_model;
results.LLH = LLH;
results.RMSE = RMSE;
results.Rsq = Rsq;
results.weights = sampleWeightAll;
results.exitflag = exitflag;
results.splits = splits;
results.ilb = ilb;
results.iub = iub;
results.lb = lb;
results.ub = ub;
results.psize = psize;

end

