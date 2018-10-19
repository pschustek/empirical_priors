classdef genEv
    %genEv Wrapper class to generate observational data for the experiments
    
    properties
    end
    
    methods(Static)
        
        %% Hierachical with binary switch between mirrored betas
        function trials = two_betas( numTrials, blockLength, bBias1, bBias2, numDist )
        %two_betas Hierarchical with binary switch between mirrored betas
        %
        %   More to follow
        
        % Sample block bias
        % Repeat same bias for all trials in block
        blockBias = repmat(binornd(1,0.5,1,numTrials/blockLength),blockLength,1);
        blockBias = blockBias(:);
        
        % Sample 'coin bias'; Stable across one block 
        pH = betarnd(bBias1*blockBias+bBias2*(1-blockBias), bBias2*blockBias+bBias1*(1-blockBias));
        
        % Sample the sample size form categorical distribution
        rep = repmat(numDist.support,numTrials,1)';
        % Categories correspond to row index
        num = rep(logical(mnrnd(1,repmat(numDist.pmf,numTrials,1)))');
        
        % Variable for coin bias
        coinBias = pH>0.5;
        % Break ties at random
        ties = pH==0.5;
        coinBias(ties) = binornd(1,0.5,sum(ties),1);
        
        % Sample observations
        meanEvidence = binornd(num,pH)./num;

        % Bundle in struct
        trials.num = [1:numTrials]';
        trials.prob = pH;
        trials.coinBias = coinBias;
        trials.meanEvidence = meanEvidence;
        trials.sampleSize = num;
        trials.blockIndex = ceil(trials.num/blockLength);
        trials.blockBias = blockBias;
        end
        
    end
end

