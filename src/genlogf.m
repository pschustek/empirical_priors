function [ fh ] = genlogf( X, w )
%genlogf Wrapper for generalized logistic sigmoid function

tmp = [X; num2cell(w')];

fh = glogf(tmp{:});

    function [ Y ] = glogf( X, varargin )
        %glogf Calculates generalized logistic function
        %   Detailed explanation goes here
        
        % Taken from:
        % http://en.wikipedia.org/wiki/Generalised_logistic_function
        
        % X is input array
        % B controls the growth rate
        % M controls lateral shift
        % A is lower asymptote
        % K is upper asymptote
        % 'nu' affects near which asymptote maximum growth occurs.
        
        % Seem to be redundant
        % Q is related to Y(0); controls the bow
        % C should be around 1
        
        % Set defaults
        optargs = {1 0 0 1 1 1 1};
        % Overwrite defaults if provided
        numvarargs = length(varargin);
        optargs(1:numvarargs) = varargin;
        
        % Place optional args in memorable variable names
        [B, M, A, K, nu, Q, C] = optargs{:};
        
        x = X(:);
        
        Y = A + (K-A)./(C+Q.*exp(-B.*(x-M))).^(1./nu);
        
        Y = reshape(Y, size(X));
    end

end
