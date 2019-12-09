%% Compliance with hierarchical task
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = 1:24;

% Load data
load('.\..\..\..\data\exp2_data.mat');

% Permutation test
nPerm = 10000;

%% Compute results
offsetDiff = nan(24,1);
perm_pval = nan(24,1);
offsetShuffled = nan(nPerm,24);
for s=subInd

    clear trials    
    trials = trialData{s};

    mEv = trials.meanEvidence;
    response = trials.confHeads;
    bBias = logical(trials.blockBias);
    
    % Offset difference
    mdl_blue = fitlm(mEv(bBias),response(bBias),'linear');
    mdl_red = fitlm(mEv(~bBias),response(~bBias),'linear');
    
    offsetDiff(s) = mdl_blue.Coefficients.Estimate(1) - mdl_red.Coefficients.Estimate(1);
    
    % Permutation test; Shuffle bBias
    res = nan(nPerm,1);
    parfor j=1:nPerm
        idx = randperm(size(response,1));
        iB = idx(1:sum(bBias));
        iR = idx(sum(bBias)+1:end);
        
        mdl_B = fitlm(mEv(iB), response(iB));
        mdl_R = fitlm(mEv(iR), response(iR));
        
        s_B = mdl_B.Coefficients.Estimate(1);
        s_R = mdl_R.Coefficients.Estimate(1);
        res(j) = s_B - s_R;
    end
    perm_pval(s) = mean(res>=offsetDiff(s));
    offsetShuffled(:,s) = res;
end

%% Save
% save('compliance_prior_task','offsetDiff','perm_pval','offsetShuffled','subInd');

%% Aggregate across participants
load('compliance_prior_task.mat'); 
[~,idx] = sort(perm_pval,'descend');

compliance = table(offsetDiff(idx), perm_pval(idx),'VariableNames',{'offset','pval'},'RowNames',num2cell(num2str(round(subInd(idx)')),2))