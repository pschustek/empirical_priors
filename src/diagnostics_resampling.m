function [ flag ] = diagnostics_resampling(w, LLH, lb, ub)
%diagnostics_resampling Heuristic checks for resampling methods on weights
%and log likelihood
%   Flags count suspicious values and should ideally be zero 

% Different from boundaries
flag.bnd = sum(sum(abs(w-ub)<1E-2 | abs(w-lb)<1E-2,2));

% Check for suspiciously large variations of weights with respect
% to overall range
flag.var = sum((quantile(w,0.9)-quantile(w,0.1))./(ub-lb) > 0.2);

% Check for outlying weights with respect to one another
ubnd = 1.5*iqr(w) + quantile(w,0.75);
lbnd = quantile(w,0.25) - 1.5*iqr(w);
flag.out = sum(sum((w > ubnd) | (w < lbnd),2));

% Check for outlying values of LLH
ubnd = 1.5*iqr(LLH) + quantile(LLH,0.75);
lbnd = quantile(LLH,0.25) - 1.5*iqr(LLH);
flag.LLH = sum(sum((LLH > ubnd) | (LLH < lbnd),2));

flag.master = sum([flag.bnd flag.var flag.out flag.LLH]);

end

