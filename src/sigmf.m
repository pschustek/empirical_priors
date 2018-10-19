function [ Y ] = sigmf( w, X )
%glogf Calculates generalized logistic function
%   Detailed explanation goes here

Y = 1./(1+exp( -(w(1) + w(2)*X + w(3)*X.^3 )));

end

