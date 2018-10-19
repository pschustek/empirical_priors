function [ out ] = rsquared( model, target )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

res = model - target;

SS_res = sum( res.^2 );
SS_tot = sum( (target-mean(target)).^2 );
out = 1-SS_res/SS_tot;

end

