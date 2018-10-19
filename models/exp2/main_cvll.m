%% Main file for calculating CV log likelihood
% For every participant and model
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\src\');

% Subselect subjects
subInd = 1:24;

% Load data
load('.\..\..\data\exp2_data.mat');

% Try to load existing structure
try 
    load prior_cvll
catch err
    if strcmp(err.identifier, 'MATLAB:load:couldNotReadFile')
       display('>>: Could not load a file')
       display('>>: Will proceed and create')
    end
end

%% Calculate
for s=subInd
    s

    %% Data
    trials = trialData{s};
    numTrials = length(trials.num);
    % block index
    block = ceil([1:numTrials]'/unique(trials.blockLength));
    % Prir belief in block tendency (variable M(b=1))
    pbEv = trials.prevConfBlockHeads;
    % Proportion of blue samples
    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    % Confidence report of the participant
    Ydat = trials.confHeads;
    
    nfolds = 5;
    
    % Create splits
    splits = [make_grouped_kfold_stratified_CV_splits( numTrials, nfolds, Ydat, block, 4587 ) ...
        make_grouped_kfold_stratified_CV_splits( numTrials, nfolds, Ydat, block, 834 ) ...
        make_grouped_kfold_stratified_CV_splits( numTrials, nfolds, Ydat, block, 237 ) ...
        make_grouped_kfold_stratified_CV_splits( numTrials, nfolds, Ydat, block, 927 ) ...
        make_grouped_kfold_stratified_CV_splits( numTrials, nfolds, Ydat, block, 04612 )];

    %% opt_opt_sig
    Xdat = 2*(trials.optConfHeads-0.5);
    f_model = @(w,X) sigmf(w, X);
    x0 = [0 6 2];
    lb = [-0.3 -1E2*[1 1]];
    ub = [0.3 1E4 1E2];
    ilb = [-0.2 -3 -3];
    iub = [0.2 12 15];
    psize = 500;
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    %opts = optimoptions('fmincon','Algorithm','interior-point','Display','iter-detailed','MaxIter',100,'FunValCheck','on','PlotFcn',@optimplotx);
    [ results ] = runCV(Xdat, Ydat, f_model, x0, lb, ub, opts, ilb, iub, psize, splits);
    cvll_prior(s).opt_opt_sig = results;
    cvll_prior(s).opt_opt_sig.flags = diagnostics_resampling(results.weights,results.LLH,results.lb,results.ub);
    
    %% opt_zmap
    f_model = @(w,X) m_map(w, X);
    Xdat = [mEv N pbEv];
    x0 = ones(1,9);
    lb = [-0.3 -15*ones(1,8)];
    ub = [0.3 15*ones(1,8)];
    ilb = [-0.2 -4*ones(1,8)];
    iub = [0.2 4*ones(1,8)];
    psize = 1000;
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    [ results ] = runCV(Xdat, Ydat, f_model, x0, lb, ub, opts, ilb, iub, psize, splits);
    cvll_prior(s).opt_zmap = results;
    cvll_prior(s).opt_zmap.flags = diagnostics_resampling(results.weights,results.LLH,results.lb,results.ub);
    
    %% tly_zmap
    f_model = @(w,X) m_map(w, X);
    Xdat = [mEv N estimate_M.tally([mEv, N, block])];
    x0 = ones(1,9);
    lb = [-0.3 -15*ones(1,8)];
    ub = [0.3 15*ones(1,8)];
    ilb = [-0.2 -4*ones(1,8)];
    iub = [0.2 4*ones(1,8)];
    psize = 1000;
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    [ results ] = runCV(Xdat, Ydat, f_model, x0, lb, ub, opts, ilb, iub, psize, splits);
    cvll_prior(s).tly_zmap = results;
    cvll_prior(s).tly_zmap.flags = diagnostics_resampling(results.weights,results.LLH,results.lb,results.ub);
    
    %% avg_zmap
    f_model = @(w,X) m_map(w, X);
    Xdat = [mEv N estimate_M.avg([mEv, N, block])];
    x0 = ones(1,9);
    lb = [-0.3 -15*ones(1,8)];
    ub = [0.3 15*ones(1,8)];
    ilb = [-0.2 -4*ones(1,8)];
    iub = [0.2 4*ones(1,8)];
    psize = 1000;
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    [ results ] = runCV(Xdat, Ydat, f_model, x0, lb, ub, opts, ilb, iub, psize, splits);
    cvll_prior(s).avg_zmap = results;
    cvll_prior(s).avg_zmap.flags = diagnostics_resampling(results.weights,results.LLH,results.lb,results.ub);
    
    save('prior_cvll', 'cvll_prior');

    %% diff_zmap
    Xdat = [mEv N block];
    f_model = @(w,X) m_map(w(2:end), [X(:,1), X(:,2), estimate_M.diff(w(1), X)]);
    x0 = [1 ones(1,9)];
    lb = [0 -0.3 -15*ones(1,8)];
    ub = [100 0.3 15*ones(1,8)];
    ilb = [1E-6 -0.2 -4*ones(1,8)];
    iub = [4 0.2 4*ones(1,8)];
    psize = 1000;
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    [ results ] = runCV(Xdat, Ydat, f_model, x0, lb, ub, opts, ilb, iub, psize, splits);
    cvll_prior(s).diff_zmap = results;
    cvll_prior(s).diff_zmap.flags = diagnostics_resampling(results.weights,results.LLH,results.lb,results.ub);
  
    %% pHT_pHT_sig
    Xdat = [mEv N block];
    f_model = @(w,X) pHT_pHT_sig(w,X);
    x0 = [5 9 0 6 2];
    lb = [0 1 -0.3 -1E2*[1 1]];
    ub = [50 60 0.3 1E4 1E2];
    ilb = [1 1 -0.2 -3 -3];
    iub = [50 60 0.2 12 15];
    psize = 1200;
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    [ results ] = runCV(Xdat, Ydat, f_model, x0, lb, ub, opts, ilb, iub, psize, splits);
    cvll_prior(s).pHT_pHT_sig = results;
    cvll_prior(s).pHT_pHT_sig.flags = diagnostics_resampling(results.weights,results.LLH,results.lb,results.ub);
    
    save('prior_cvll', 'cvll_prior');
    
    % Display output
    msg = sprintf('>>: Iteration %d completed', s);
    disp(msg);

end