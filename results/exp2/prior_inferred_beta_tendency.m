%% Beta distribution of block tendency fitted to data
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\src\');

% Subselect subjects
subInd = 1:24;

% Load data
load('.\..\..\models\exp2\prior_ml.mat')

%% 
pHT_pHT_weights = nan(numel(subInd),5);
j = 0;
for s=subInd
    j = j + 1;
    pHT_pHT_weights(j,:) = ml_prior(s).pHT_pHT_sig.weights;
end

%% pHT_pHT_sig prior fit
t0 = pHT_pHT_weights(:,2);
h0 = pHT_pHT_weights(:,1) + t0;

muExp = h0./(h0+t0);
varExp = h0.*t0./((h0+t0).^2.*(h0+t0+1));

muMod = 14/(14+9);
varMod = 14.*9./((14+9).^2.*(14+9+1));

fprintf('- [result] median = %.4f\n', median(muExp));

% Different from optimal
[p,~,stats] = signrank(muExp,muMod,'tail','left');
fprintf('- [result] p-value = %.4e\n', p);