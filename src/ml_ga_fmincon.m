function [ LLH, Rsq, fweights, exitflag, g_llh, w ] = ml_ga_fmincon( Xdat, Ydat, func, x0, lb, ub, opts, ilb, iub, psize )
%ml_ga_fmincon Find ML solution for given model
%   Detailed explanation goes here

% Convert to minimization problem
% Objective funtion computes log likelihood of the responses
fOpt = @(wFit)  -f_obj_adapt( func(wFit, Xdat), Ydat );

% Genetic algorithm
% rng_settings = rng('shuffle','twister');
ipb = [ilb; iub];

% Setting for genetic algorithm
opts_ga = gaoptimset('UseParallel','always','PlotFcn',{@gaplotbestf,@gaplotbestindiv,@gaplotdistance,@gaplotscorediversity},'Display','iter',...
    'PopInitRange',ipb,'MigrationDirection','both','PopulationSize',psize,'TimeLimit',700,'StallGenLimit',4,'StallTimeLimit',800,...
    'TolFun',5E-5,'Generations',40);
[w, fval] = ga(fOpt,length(x0),[],[],[],[],lb,ub,[],opts_ga);

% Log likelihood
g_llh = -fval;

% Refine ML estimate with regression
[fweights,fval,exitflag] = fmincon(fOpt,w,[],[],[],[],lb,ub,[],opts);

Ymodel = func(fweights, Xdat);

LLH = -fval; 
Rsq = rsquared( Ymodel, Ydat );

end

