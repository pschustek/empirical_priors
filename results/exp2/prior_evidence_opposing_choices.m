%% Measure whether participants go against momentary evidence 
clear variables
close all

%% Prepare
addpath('.\..\..\src\')
load('.\..\..\data\exp2_data.mat')

%% Prepare data
subInd = [1:20 22:24];

rng('shuffle');

%% 
numOpp = nan(max(subInd),1);
numOppOpt = nan(max(subInd),1);
numOppmEv = nan(max(subInd),1);
for s=subInd

    clear trials    
    trials = trialData{s};

    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    blockLength = trials.blockLength;
    numTrials = size(trials.num,1);
    blockBias = trials.blockBias;

    pbConf = trials.prevConfBlockHeads;
    response = trials.confHeads;
    optRes = trials.optConfHeads;
    % Augment with different noise instances (null hypothesis)
    fac = 10000;
    nullHyp = trunc_gauss_sample(repmat(mEv,fac,1), 0.1, 0, 1);

    % Align with real block bias (symmetric pattern)
    block = 2*(trials.blockBias-0.5);
    response = block.*(response-0.5) + 0.5;
    optRes = block.*(optRes-0.5) + 0.5;
    mEv = block.*(mEv-0.5) + 0.5;
    pbConf = block.*(pbConf-0.5) + 0.5;
    nullHyp = repmat(block,fac,1).*(nullHyp-0.5) + 0.5;
    
    % Sample evidence opposing real block tendency; Main contribution comes
    % from those trials
    mask = mEv < 0.5;
    % Number of choices against momentary evidence
    numOpp(s) =  sum(response(mask)>0.5)/numTrials;
    numOppOpt(s) =  sum(optRes(mask)>0.5)/numTrials;   
    numOppmEv(s) =  sum(nullHyp(repmat(mask,fac,1))>0.5)/(numTrials*fac);
    
end

%% Aggregate across participants
% Reference: betacdf(0.5,14,9) = 0.1431
nExp = numOpp(subInd);
nOpt = numOppOpt(subInd);
nNull = numOppmEv(subInd);

[pMod,~,statsMod] = signrank(nExp,nOpt,'tail','left');
fprintf('- [result] Fewer than optimal, p-value = %.4f\n', pMod);

[pNull,~,statsNull] = signrank(nExp,nNull,'tail','right');
fprintf('- [result] More than noise, p-value = %.4f\n', pNull);