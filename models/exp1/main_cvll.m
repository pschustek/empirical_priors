%% Main script for calculating cross-validation log likelihood (CVLL) for Experiment 1
% For every participant and model

%% Setup
clear variables
close all

% Add path to auxiliary functions
addpath('.\..\..\src\');

% Load data
load('.\..\..\data\exp1_data.mat');

% starting points and boundaries for distorsion parameters
x0_sig = [0 1 0];
lb_sig = [-0.3 -1E2*[1 1]];
    ub_sig = [0.3 1E2*[1 1]];
    ilb_sig = [-0.2 -3 -3];
    iub_sig = [0.2 12 15];


sessionIndices = [1:24];

% Try to load existing structure
try
    load basic_cvll
catch err
    % Just to let us know; Execution continues
    if strcmp(err.identifier, 'MATLAB:load:couldNotReadFile')
        display('>>: Could not load a file')
        display('>>: Will proceed and create')
    end
end

%% Calculate CVLL for every participant
for s=sessionIndices
    s
    
    %% Data
    clear trials
    trials = trialData{s};
    numTrials = length(trials.num);
    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    NB = mEv .* N;
    NR = N-NB;
    % Confidence report of the ideal observer
    optConfH = trials.optConfHeads;
    % Confidence report of the participant
    Ydat = trials.confHeads;
    
    nfolds = 5;
    % Generate CV folds with fixed seed of random number generator
    splits = [makeCVsplits_kfold_strat( numTrials, nfolds, Ydat, 621 ) ...
        makeCVsplits_kfold_strat( numTrials, nfolds, Ydat, 358 ) ...
        makeCVsplits_kfold_strat( numTrials, nfolds, Ydat, 4673 ) ...
        makeCVsplits_kfold_strat( numTrials, nfolds, Ydat, 32 ) ...
        makeCVsplits_kfold_strat( numTrials, nfolds, Ydat, 9276 )];
    
    % Settings of the genetic algorithm hard-coded in respective file
    % ml_ga_fmincon()
    
    %% Mean evidence model
    Xdat = 2*(mEv-0.5);     % standardize input
    f_model = @(w,X) sigmf( w, X );
    % Setting for fmincon
    x0 = x0_sig; [0 6 2];
    lb = lb_sig; [-0.3 -1E2*[1 1]];
    ub = ub_sig; [0.3 1E2 1E2];
    ilb = ilb_sig; [-0.2 -3 -3];
    iub = iub_sig; [0.2 12 15];
    psize = 400;
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    [ results ] = runCV(Xdat, Ydat, f_model, x0, lb, ub, opts, ilb, iub, psize, splits);
    cvll_basic(s).mEv = results;
    cvll_basic(s).mEv.flags = diagnostics_resampling(results.weights,results.LLH,results.lb,results.ub);
    
    %% Difference model
    nH = mEv.*N;
    D = nH - (N-nH);
    f_model = @(w,X) sigmf( w, X );
scalef = max(N).^[0 -1 -3]; % scaling factor

   x0 = x0_sig .*scalef;
    lb = lb_sig .*scalef; 
    ub = ub_sig .*scalef; 
    ilb = ilb_sig .*scalef;
    iub = iub_sig .*scalef; 
    psize = 400;
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    [ results ] = runCV(D, Ydat, f_model, x0, lb, ub, opts, ilb, iub, psize, splits);
    cvll_basic(s).diff = results;
    cvll_basic(s).diff.flags = diagnostics_resampling(results.weights,results.LLH,results.lb,results.ub);
    
   
    
    %% Ideal observer model
    Xdat = 2*(optConfH-0.5);
    f_model = @(w,X) sigmf( w, X );
  x0 = x0_sig;
    lb = lb_sig; 
    ub = ub_sig; 
    ilb = ilb_sig;
    iub = iub_sig; 
    psize = 400;
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    [ results ] = runCV(Xdat, Ydat, f_model, x0, lb, ub, opts, ilb, iub, psize, splits);
    cvll_basic(s).opt = results;
    cvll_basic(s).opt.flags = diagnostics_resampling(results.weights,results.LLH,results.lb,results.ub);
    
       
    %%
    save('basic_cvll', 'cvll_basic');
    
    % Display output
    msg = sprintf('>>: Iteration %d completed', s);
    disp(msg);
    
end