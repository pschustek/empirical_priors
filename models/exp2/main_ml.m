%% Main file for estimating parameters by ML
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
    load prior_ml
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
    block = ceil([1:numTrials]'/unique(trials.blockLength));
    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    Ydat = trials.confHeads;
    
    %% pHT_pHT_sig
    Xdat = [mEv N block];
    f_model = @(w,X) pHT_pHT_sig(w,X);
    x0 = [5 9 0 6 2];
    lb = [0 1 -0.3 -1E2*[1 1]];
    ub = [50 60 0.3 1E4 1E2];
    ilb = [1 1 -0.2 2 -3];
    iub = [50 60 0.2 12 15];
    psize = 1200;
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    [ LLH, Rsq, fweights, exitflag, g_llh, g_w ] = ml_ga_fmincon( Xdat, Ydat, f_model, x0, lb, ub, opts, ilb, iub, psize );
    ml_prior(s).pHT_pHT_sig.LLH = LLH;
    ml_prior(s).pHT_pHT_sig.Rsq = Rsq;
    ml_prior(s).pHT_pHT_sig.weights = fweights;
    ml_prior(s).pHT_pHT_sig.exitflag = exitflag;
    ml_prior(s).pHT_pHT_sig.g_llh = g_llh;
    ml_prior(s).pHT_pHT_sig.g_w = g_w;
    ml_prior(s).pHT_pHT_sig.lb = lb;
    ml_prior(s).pHT_pHT_sig.ub = ub;
    ml_prior(s).pHT_pHT_sig.ilb = ilb;
    ml_prior(s).pHT_pHT_sig.iub = iub;
    ml_prior(s).pHT_pHT_sig.psize = psize;

    save('prior_ml', 'ml_prior');

    % Display output
    msg = sprintf('>>: Iteration %d completed', s);
    disp(msg);

end
%% Save results
save('prior_ml', 'ml_prior'); 