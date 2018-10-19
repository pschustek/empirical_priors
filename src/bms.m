function [ a ] = bms( num, M )
%bms Bayesian model selection
%   Implemented after 'Bayesian model selection for group studies', Stephan
%   et al., NeuroImage, 2009


% 1st step
alph0 = ones(1,size(M,2));
a = alph0;

% Iteration
for i=1:num
    u = nan(size(M));
    denom = nan(size(M,1),1);
    for n=1:size(M,1)       % loop over participants
        u(n,:) = exp( M(n,:) + psi(a) - psi(sum(a)) );

        denom(n) = sum(u(n,:));

    end

    b = nan(1,size(M,2));
    for k=1:size(M,2)       % loop over models
        b(k) = sum(u(:,k)./denom);
    end

    a = alph0 + b;
end

end

