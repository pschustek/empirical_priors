clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\src\');

% Subselect subjects
subInd = [1:20 22:24];

% Load data
load('.\..\..\models\exp2\model_evidence')

M = M(subInd,:);
rowNames = {rowNames{subInd}};

%tM = array2table(M,'RowNames',rowNames,'VariableNames',fnames)

%% Select (order of) models 
list = {'pHT_pHT_sig','opt_opt_sig'};
modelLabel = {'pHT_sig','opt_sig'};

IDX = nan(length(list),1);
for m=1:length(list)  
    IDX(m) = find(strcmp(fnames, list{m}));
end

%% Pairwise BMS comparison
numIter = 50;
numModels = length(IDX);

medDiff = nan(1,numModels-1);
exceed = nan(1,numModels-1);
expected = nan(1,numModels-1);
expectedCI = nan(2,numModels-1);

for j=2:numModels
    comp = [IDX(1) IDX(j)];         % compare only pairs
    alpha = bms(numIter,M(:,comp));
    
    fh = @(q) bms(numIter,q);
    
    % Bootstrapped BMS
    bs = bootstrp(1000,fh,M(:,comp));
    R = bs./repmat(sum(bs,2),1,2);
    
    % Expected likelihood of model k
    r = alpha./sum(alpha);
    expected(1,j-1) = r(1);
    expectedCI(:,j-1) = quantile(R(:,1),[0.05 0.95])';
    
    % Standard difference
    medDiff(1,j-1) = -median(dHart(diff(M(:,comp),1,2)));
    
    % Beta-marginal
    exceed(1,j-1) = betacdf(0.5,alpha(1),sum(alpha)-alpha(1),'upper');       % actual exceedance probability: 1 vs. the other
end

%% Result
fprintf('- [result] exceedance probability = %.6f\n', exceed(1));
