function [ LL ] = f_obj( Ymodel, Y, varargin )
%f_obj Model of the response noise
%   Use of normal distribution

% Width of constant Gaussian noise
if nargin == 3
    RMSD = varargin{1};
else
    % Estimate noise
    RMSD = sqrt(mean( (Ymodel-Y).^2 ));
end

eps = 1.34E-4;
LLP = (1-eps)*normpdf(Y,Ymodel,RMSD) + eps;
LL = sum(log(LLP));

assert(isnan(LL)==0);
assert(isreal(LL)==1);

end

