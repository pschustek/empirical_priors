function Y = m_map(w, X)
%m_map Flexible map of M(b) and sample statistics onto the response
%   Only wrapper for nested function
%

tmp = [{X(:,1); X(:,2); X(:,3)}; num2cell(w')];
Y = m_map_calc(tmp{:});

    function [ Y ] = m_map_calc( mEv, N, pbEv, varargin)
        %m_map Summary of this function goes here
        %   Detailed explanation goes here
        %
        
        % Set defaults
        optargs = num2cell(zeros(1,9));
        % Overwrite defaults if provided
        numvarargs = length(varargin);
        optargs(1:numvarargs) = varargin;
        % Place optional args in memorable variable names
        [q1, q2, q3, q4, q5, q6, q7, q8, q9] = optargs{:};
        
        % Normalize
        pbEv = 2*(pbEv - 0.5);
        mEv = 2*(mEv - 0.5);
        
        % Standardize input values (mEv and pbEv are already centered)
        Z = q1 + q2*mEv/0.49 + q3*mEv.*N/3.35 + q4*pbEv/0.57 + q5*mEv.^3/0.34 + ...
            q6*mEv.^3.*N/1.98 + q7*pbEv.*N/4.24 + q8*pbEv.^3/0.41 + q9*pbEv.^3.*N/3.02;
        
        Y = 1./(1+exp(-Z));
    end

end
