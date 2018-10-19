function [out] = trunc_gauss_pdf(q, mu, sig, lb, ub)
%trunc_gauss_pdf Truncated normal distribution
%   Detailed explanation goes here

out = normpdf(q,mu,sig)./(normcdf(ub,mu,sig)-normcdf(lb,mu,sig));

end

