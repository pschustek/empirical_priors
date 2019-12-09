function [ fweights, R2 ] = fit_prior_weights_fixed( Xdat, Ydat, numWeights )
%fit_prior_weights Influence of each previous trial on the current response
%   For different number of previous trials

num = unique(numWeights);

%% Settings
S = 2;
L = -10;
U = 10;
opts = optimoptions('fmincon','Algorithm','active-set','Display','none');
%opts = optimoptions('fmincon','Algorithm','interior-point','Display','iter-detailed','MaxIter',100,'FunValCheck','on','PlotFcn',@optimplotx);

%% Regression model
% Separately for different inBlock trials (different number of predictors)
ngroups = numel(num);
nreg = size(Xdat,2); % number of regressors
R2 = nan(ngroups,1);
fweights = nan(ngroups,max(num));
k = 0;

% remove 5 first trials
Xdat(1:ngroups,:) = [];
Ydat(1:ngroups) = [];
numWeights(1:ngroups) = [];

for j=unique(numWeights)'
    k = k + 1;
    X = Xdat(numWeights==j,:);
    Y = Ydat(numWeights==j);

    f_min = @(w) -f_obj( prior_trialWeights(w,X),Y);
    [w,~,~] = fmincon(f_min,ones(1,nreg)*S,[],[],[],[],ones(1,nreg)*L,ones(1,nreg)*U,[],opts); 
    R2(k) = rsquared( prior_trialWeights( w, X ), Y );
    fweights(k,:) = w;
end

end

function [ output ] = prior_trialWeights( w, X )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

reg = sum(repmat(w,size(X,1),1).*X,2);

output = 1./(1+exp( -(reg)));

end

