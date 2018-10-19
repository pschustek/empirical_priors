classdef opt_inf
    %opt_inf Various infererences for different inputs
    %   Assumes the two-betas hierarchical structure with binary block bias switch.
    
    properties
    end
    
    methods(Static)
        %% Complete inferernce given sample data
        function [ confHeads, bEv, pbEv, confBlockHeads ] = all( nH, N, block, bBias1, bBias2 )
        %all Inference of confidence and messages of all trials in a block
        %   confHeads = Confidence in airplane majority for each trial
        %   bEv = Evidence for block variable b; message m(b=1)
        %   pbEv = Prior belief in the contextual variable b; M(b=1)
        %   confBlockHeads = Belief in the value of b after observing the sample
        
        numTrials = length(nH);
        nT = N - nH;
        
        % Fast computation with betacdf based on beta-function (normalization factor of beta distribution)
        bEv = 1./( 1 + gamma(nH+bBias2).*gamma(nT+bBias1)./(gamma(nH+bBias1).*gamma(nT+bBias2)) );

        % Calculate confidence in block bias (accumulative)
        un = unique(block)';
        mask = repmat(block,1,length(un))==repmat(un,numTrials,1);
        tmp = repmat(bEv,1,size(mask,2));
        blockLength = unique(sum(mask));
        assert(length(blockLength)==1);
        bEv_H_inBlock = reshape(tmp(mask),blockLength,[]);          
        bEv_T_inBlock = reshape(1-tmp(mask),blockLength,[]);
        % Column-wise cumulative product
        cumHev = cumprod(bEv_H_inBlock,1);
        cumTev = cumprod(bEv_T_inBlock,1);
        % Renormalization
        normCol = cumHev + cumTev;
        confBlockHeads = reshape(cumHev./normCol, numTrials,1);

        % Confidence in block bias on previous trial
        mask = false(numTrials,1);
        mask(1:blockLength:end) = 1;
        pbEv = [nan; confBlockHeads(1:end-1)];
        pbEv(mask) = 0.5;

        %% Calculate confidence in H trial bias (momentary)
        % Calculate contribution from the two possible block bias
        % structures with the correct normalization scheme
        confHeads = ( pbEv.*betacdf(0.5,nH+bBias1,nT+bBias2,'upper').*gamma(nH+bBias1).*gamma(nT+bBias2) + ...
            (1-pbEv).*betacdf(0.5,nH+bBias2,nT+bBias1,'upper').*gamma(nH+bBias2).*gamma(nT+bBias1) ) ./ ...
            ( pbEv.*gamma(nH+bBias1).*gamma(nT+bBias2) + (1-pbEv).*gamma(nH+bBias2).*gamma(nT+bBias1) );
        
        end
        
        %% Inference 
        function [ confHeads ] = confH_bEv( nH, N, pbEv, bBias1, bBias2 )
        %confH_bEv Confidence in 'H' bias given pbEv
        %   Given momentary observations and a prior belief in the block
        %   tendency pbEv

        % Make column vectors
        nH = nH(:);
        N = N(:);
        
        pbEv = pbEv(:);

        % Pick maximum size
        numTrials = max([length(nH) length(N) length(pbEv)]);

        % Scalar expansion
        nH = nH.*ones(numTrials,1);
        N = N.*ones(numTrials,1);
        nT = N - nH;
        pbEv = pbEv.*ones(numTrials,1);
        
        % Calculate contribution from the two possible block bias
        % structures with the correct normalization scheme
        confHeads = ( pbEv.*betacdf(0.5,nH+bBias1,nT+bBias2,'upper').*gamma(nH+bBias1).*gamma(nT+bBias2) + ...
            (1-pbEv).*betacdf(0.5,nH+bBias2,nT+bBias1,'upper').*gamma(nH+bBias2).*gamma(nT+bBias1) ) ./ ...
            ( pbEv.*gamma(nH+bBias1).*gamma(nT+bBias2) + (1-pbEv).*gamma(nH+bBias2).*gamma(nT+bBias1) );
        
        end
        
        %% Momentary evidence for block bias
        function [ bEv ] = bEv( nH, N, bBias1, bBias2 )
        %bEv Momentary evidence for the block bias from each trial
        %   Independent of prior or previous trials
        %   Just marginalizes out the parameter 'u' representing the within-trial
        %   majority.

        nT = N - nH;
        
        % Fast computation with betacdf based on beta-function (normalization factor of beta distribution)
        bEv = 1./( 1 + gamma(nH+bBias2).*gamma(nT+bBias1)./(gamma(nH+bBias1).*gamma(nT+bBias2)) );

        end
        
        %% Basic task: Confidence in H as a function of the prior (fast)
        function [ confHeads ] = basic_confH( nH, N, a0, b0 )
            %basic_confH Confidence in blue airplane majority of Experiment 1
            %   Detailed explanation goes here

            nT = N - nH;
            
            confHeads = 1 - betacdf(0.5, nH+a0, nT+b0);
        end
        

    end
    
end

