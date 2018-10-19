function output = pHT_pHT_sig(w, X)
%UNTITLED pHT_pHT_sig Probabilistic inference with different block tendency
%   Fitting the parameters of the Beta distribution modeling the block
%   tendency.
%   Furthermore, distortions on the estimate (when reporting) are accounted
%   for with a logistic sogmoidal mapping on the output.

confidence = opt_inf.all(X(:,1).*X(:,2), X(:,2), X(:,3), sum(w(1:2)), w(2));
assert(sum(isnan(confidence))==0);   % function fails for high values

% Standardize
Y = 2*(confidence-0.5);

% Logistic sigmoid
output = sigmf( w(3:5), Y ); 
end

