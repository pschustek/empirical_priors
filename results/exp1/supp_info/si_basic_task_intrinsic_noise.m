%% Estimate behavioral variance for fixed input 
% For fixed sample input

clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = [1:24];

% Load data
load('.\..\..\..\data\exp1_data.mat')

%% Compute standard deviation
var_intrinsic.sub = nan(numel(subInd),1);
var_intrinsic.num = nan(numel(subInd),1);

k = 0;
for s=subInd

    k = k + 1;
    
    clear trials    
    trials = trialData{s};

    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    response = trials.confHeads;

    setN = unique(N);
    
    % Assume symmetry
    % Collapse onto distance from decision boundary; 'blue' samples
    sg = 2*((mEv >= 0.5)-0.5);                  % sign
    mQv = abs(mEv-0.5) + 0.5;   
    % Mirror responses with respect to decision boundary for 'red' samples
    response = sg.*(response-0.5) + 0.5;        
    
    sub = [];
    for n=setN'
        mask = n==N;
        pattern = ([ceil(n/2):n]/n)';       % pool 'red' and 'blue' samples of same fraction
        mask = mQv==pattern' & mask;
        num = sum(mask);
        check = mask(:,num>=10);
        for j=1:size(check,2)
            D = response(check(:,j));
            sub = [sub; (D-mean(D)).^2];
        end
    end
    
    var_intrinsic.sub(k) = mean(sub);
    var_intrinsic.num(k) = numel(sub);
end


%% Compute sigma for each participant
sigm = sqrt(var_intrinsic.sub);
bs = bootstrp(2000,@median,sigm);

fprintf('- [result] estimated SD = %.4f, 95%%-CI = (%.4f, %.4f)\n', median(sigm), quantile(bs,[0.05 0.95]));