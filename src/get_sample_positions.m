function [gridX, gridY, color] = get_sample_positions(nH, N, varargin)

% Set random number generator
if nargin >=3
    rng(varargin{1});
else
    rng('shuffle');
end

% Array of individual coin tosses
input = zeros(1,N);
input(randsample(N,nH)) = 1;     % without replacement
% Hard-coded parameter
maxN = 13;

numCoins = size(input,2);

% Y-shuffled
extra = 3 + randi(4,1);
maxN = maxN + extra;       % add random component to number of grid positions

% Two slightly seperated groups
margin = 0.05;
xGrid = linspace(margin,1-margin,maxN);
yGrid = linspace(margin,1-margin,maxN);                    % spread along y

% Sample (minimal) bounds for the two groups
% They contain an amount equal of their number plus a margin of 2
input = sort(input,'descend');
% First grid-point + length + buffer
lb = 1 + sum(input) + 1;
% Last grid-point - length - buffer
ub = maxN - (sum(~input)+1);
% The resulting separator is randomly split
sep = datasample(lb:ub,1);
% The margin are re-established (-2, +2) around the separator
% The grid positions are then subsampled for each group
gIdx = [datasample(1:(sep-2),sum(input),'Replace',false) datasample((sep+2):maxN,sum(~input),'Replace',false)];

% 'Blue' is either left or right, by flipping the direction of the xGrid
if binornd(1,0.5)
    xGrid = -(xGrid-0.5) + 0.5;
end
gridX = xGrid(gIdx);
gridY = yGrid(randsample(maxN,numCoins,false));

% Center
gridX = gridX - 0.5;
gridY = gridY - 0.5;

color = [1-input; zeros(1,numCoins); input];

end