function [ Y ] = sigmf2( w, X )
%glogf Calculates generalized logistic function from inverse logistic
%   Detailed explanation goes here

LLR = log(X./(1-X));
Y = 1./(1+exp( -(w(1) + w(2)*LLR + w(3)*LLR.^3 )));

end

