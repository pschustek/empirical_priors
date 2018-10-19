function [ confHeads, mEv, N, blockBias, pbB, bEv, coinBias ] = sim_opt_inf( numSim, varargin  )
%sim_opt_inf Summary of this function goes here
%   Detailed explanation goes here

% Set defaults
optargs = num2cell([14 9 5]);
% Overwrite defaults if provided
numvarargs = length(varargin);
optargs(1:numvarargs) = varargin;
% Place optional args in memorable variable names
[bBias1, bBias2, blockLength] = optargs{:};

%% Setup
% Distribution over sample size
numDist.support = 3:11;
numMap = @(n) 0*n + 1;
numDist.pmf = numMap(numDist.support)/sum(numMap(numDist.support));

% Split all trials into chunks
chunk = 1000;
numIt = ceil(numSim/chunk);

confHeads = nan(numIt*chunk,1);
bEv = nan(numIt*chunk,1);
pbB = nan(numIt*chunk,1);
mEv = nan(numIt*chunk,1);
blockBias = nan(numIt*chunk,1);
coinBias = nan(numIt*chunk,1);
N = nan(numIt*chunk,1);

for j=1:numIt
    % Indices of chunk
    idx = [(j-1)*chunk+1:j*chunk];
    
    % Make data
    expSetup = genEv.two_betas( chunk, blockLength, bBias1, bBias2, numDist );

    num = expSetup.sampleSize; H = expSetup.meanEvidence.*num; block = expSetup.blockIndex;
    mEv(idx) = H./num;
    blockBias(idx) = expSetup.blockBias;
    N(idx) = num;
    coinBias(idx) = expSetup.coinBias;

    % Get model predictions
    [ confHeads(idx), bEv(idx), pbB(idx)] = opt_inf.all( H, num, block, bBias1, bBias2 );
end

end