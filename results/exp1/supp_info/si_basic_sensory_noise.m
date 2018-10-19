clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = [1:24];

% Load data
load('.\..\..\..\data\exp1_data.mat');

%% 
for s=subInd
    clear trials
    trials = trialData{s};
    
    % Psychometric curve by fitting cumulative normal distribution
    ferr = @(q) mean((normcdf(trials.meanEvidence, 0.5, q) - trials.decision).^2);
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    w = fmincon(ferr,0.1,[],[],[],[],0,10,[],opts);
    sub(s).ncdf_std = w;
    
%     clf;
%     hold on
%     scatter(trials.meanEvidence, trials.decision);
%     q = 0:0.01:1;
%     plot(q,normcdf(q,w(1),w(2)), 'r');
end

%% Aggregate across participants
Q = [sub.ncdf_std]';
sdMedian = median(Q);
bs = bootstrp(1000,@median,Q);
quantile(bs,[0.05 0.95]);

fprintf('- [result] estimated Gaussian SD = %.4f, 95%%-CI = (%.4f, %.4f)\n', median(Q), quantile(bs,[0.05 0.95]));