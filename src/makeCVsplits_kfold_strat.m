function [ splits, sett ] = makeCVsplits_kfold_strat( numTrials, kfold, Y, varargin )
%makeCVsplits_kfold_strat Split data for CV while maintaining an
%approximately balanced statistics of the target variable
%   Detailed explanation goes here

% Use seed if provided
if nargin == 4
    seed = varargin{1}; 
    rng(seed);
    sett = rng;
else
    % Random seed otherwise
    rng('shuffle');
    sett = rng;
end

% Integer divisions of the number of trials
assert(mod(numTrials,kfold)==0);

%% Make splits
% Divide dependent variable into numFold bins
% nstrat and kfold must both divide numTrials
% %%%
% Use factor(numTrials) to find integer partition
cp = [cumprod(factor(numTrials/kfold)) factor(numTrials/kfold)];
cp = sort(unique(cp));
% Divisors close to 8 are preferred that satisfy the requirements
indices = find((cp>=2) & (cp<cp(end)) & (cp<=11));
% preference to be close to 8
distance = abs(cp(indices) - 8);
% Find minimum
indexOfMin = indices(min(distance)==distance);
% Pick the closest divisor
nstrat = cp(indexOfMin);

assert(mod(numTrials,nstrat)==0);       % Consistency check
numInStrat = numTrials/nstrat;
beg = 1+([1:nstrat]-1)*numInStrat;      % Starting indeces of slices
eend = [1:nstrat]*numInStrat;           % End indeces of slices
[~, idx] = sort(Y);                     % With respect to the sorted responses

% Make random but approx. balanced slices of data (one from each stratum)
S = nan(numTrials/nstrat,nstrat);
for j=1:nstrat
    % Slices indexed by column
    D = idx(beg(j):eend(j));
    % Permute trials within each slice
    I = randperm(numInStrat)';
    S(:,j) = D(I);
end

% Split S vertically according to folds
sz = size(S,1)/kfold;
% Test whether integer division is fulfilled
assert(mod(sz,round(sz))==0)    
beg = 1+([1:kfold]-1)*sz;
eend = [1:kfold]*sz;
% Logical array to be filled; folds in columns
splits = nan(numTrials,kfold);
for j=1:kfold
    % Aggregate across slices
    % The j-th split (out of k) of both slices   
    tmp = S(beg(j):eend(j),:);      
    % Find the logical indices
    % These trials will be the test set for the j-th split
    mask = any(repmat([1:numTrials]',1,numel(tmp))==repmat(tmp(:)',numTrials,1),2);
    splits(mask,j) = 0;         % 0 stands for test set trial
    splits(~mask,j) = 1;        % 1 stands for training set trial
end
splits = logical(splits);

% Check for correct number of training trials in each split
assert(sum(sum(splits)==(numTrials*(kfold-1)/kfold))==kfold);
% Check for (kfold-1):kfold use of each fold for training across folds
assert(sum(~(sum(splits,2)==(kfold-1)))==0);

% Check for approximate balance of mean()
tmp = repmat(Y,1,kfold);
m1 = mean(reshape(tmp(splits),numTrials*(kfold-1)/kfold,[]));
m2 = mean(reshape(tmp(~splits),numTrials/kfold,[]));
% Mean relative absolute error of the response mean
mq = mean(abs(m1-m2)./m1);

% Check for approximate balance of std()
s1 = std(reshape(tmp(splits),numTrials*(kfold-1)/kfold,[]));
s2 = std(reshape(tmp(~splits),numTrials/kfold,[]));
% Mean relative absolute error of the response SD
sq = mean(abs(s1-s2)./s1);

end

