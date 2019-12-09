function [ LL, S ] = f_obj_adapt( Ymodel, Y, varargin )
%f_obj_adapt Model of the response noise
%   Use of truncated normal distribution

eps = 1.34E-4;  realmin; % normpdf(4,0,1)

if nargin>=3
    % Use passed value as function argument
    assert(length(varargin)==1)
    S = varargin{1};
end
    
% Estimate S
if nargin==2
    res = Y-Ymodel;
    
    % Remove outliers from calculation of S
    % With respect to residuals
    lb = quantile(res,0.25) - 3*iqr(res);
    ub = quantile(res,0.75) + 3*iqr(res);
    trimmed = res((res>=lb)&(res<=ub));
    S = sqrt(mean(trimmed.^2));
    S = max(S, 0.01); % minimal value for noise to avoid singular solutions explaining boundary responses
end

% Calculate LL for all values (even trimmed)
LLP = (1-eps)*trunc_gauss_pdf(Y,Ymodel,S,0,1) + eps;

% Remove nan values
LLP(isnan(LLP)) = [];
% Log likelihood of all considered points
LL = sum(log(LLP));


end

