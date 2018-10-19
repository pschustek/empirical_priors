function [ out ] = dHart( LLR )
%dHart Convert log probability ratio to decibans (or Hartleys)
%   log probability ratio as input

out = log10(exp(LLR))*10;

end

