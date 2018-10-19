function [ splits, sett ] = make_grouped_kfold_stratified_CV_splits( numTrials, kfold, Y, groups, varargin )
%make_grouped_kfold_stratified_CV_splits Respect grouping when making CV
%splits
%   Detailed explanation goes here

if nargin == 5
    seed = varargin{1}; 
    rng(seed);
    sett = rng;
else
    rng('shuffle');
    sett = rng;
end

numGroups = numel(unique(groups));
% Maximum amount of integer divisible indices (groups)
nmod = mod(numGroups,kfold);

% Consistency check for constant group size
assert(mod(numTrials,numGroups)==0);

%% Split only groups of trials
% Logical trial indices for each group (columns)
group2Idx = nan(numTrials,numGroups);   
% Mean response for each group of trials
meanYgrp = nan(numGroups,1);
for b=1:numGroups
    group2Idx(:,b) = groups == b;
    meanYgrp(b) = mean(Y(groups==b));
end

% Variable used for stratification; Symmetry
Yq = abs(meanYgrp-0.5)+0.5;
% Select group of remaining (modulo) groups
R = randi(numGroups,1,nmod);
% Mask for all other groups
mask = ~any([1:numGroups]'==R,2);
% Use method on group indices
% Pass current state of random number generator 
x = makeCVsplits_kfold_strat( numGroups-nmod, kfold, Yq(mask), rng() );

% Merge remaining modulo groups
S = nan(numGroups,kfold);
% Fill with already partitioned groups
S(mask,:) = x;
% Proportional fraction of test trials
% Number of training/test trials for each fold should be the same
I = round(nmod/kfold);
S(R(1:I),:) = 0;                % test
S(R(min(I+1,nmod):end),:) = 1;  % training

S = logical(S);

%% Make splits
% Map from groups to trials
splits = zeros(numTrials,kfold);
for i=1:kfold
    splits(:,i) = any(group2Idx(:,S(:,i)),2);
end
splits = logical(splits);

% Number of training trials is an integer multiple of (constant) group size for each fold 
assert(prod(mod(sum(splits),numTrials/numGroups))==0)

end

