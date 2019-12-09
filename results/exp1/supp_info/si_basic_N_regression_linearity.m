% Supp Fig. 3c
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = 1:24;

% Load data
load('.\..\..\..\data\exp1_data.mat');

% Load model fits
load('.\..\..\..\models\exp1\basic_cvll.mat');

modelnames = {'mEv','diff','opt'}; % models
model_labels = {'ratio', 'difference','optimal'};
n_models = length(modelnames);

nPerm = 10000;

%% Compute pattern
% Takes ~ 10 min
tic;
for s=subInd % for each subject
    
    clear trials
    trials = trialData{s};
    fprintf('- computing pattern for participant %d\n', s);
    
    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    numTrials = size(trials.num,1);
    
    response = trials.confHeads; % subject response
    optRes = trials.optConfHeads; % optimal model response
    
    %% Fit N
    setN = unique(N); % sample sizes used
    x0 = [0 0];
    lb = -1E2*[1 1];
    ub = 1E2*[1 1];
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    
    % Normalize mEv
    mEv = 2*(mEv-0.5);
    
    
    %  for j=1:length(setN) % for each sample size
    %      mask = setN(j)==N; % select trials with given sample size
    %    ff = @(w,X) genlogf(X,w); % logistic function (only parameter is slope)
    ff = @(w,X, N) 1./(1+exp(-w(1)*N.^w(2).*X));
    
    % Participant
    fObj = @(w) -f_obj( ff(w,mEv, N), response );
    [sub(s).w,~,exitflag] = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
    
    % Model
    fObj = @(w) -f_obj( ff(w,mEv, N), optRes );
    [sub(s).wOpt,~,exitflag] = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
    
    % slope for fitted models
    for m=setdiff(1:n_models,5)
        strk =  cvll_basic(s).(modelnames{m}); % model
        fObj = @(w) -f_obj( ff(w,mEv, N), strk.Y );
        [sub(s).wMod(m,:),~,exitflag] = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
        
    end
    
    %         clf;
    %         hold on
    %         scatter(mEv(mask),response(mask));
    %         q = -1:0.05:1;
    %         plot(q,ff(slope(j),q),'r');
    %         LL = logLikeGauss( ff(x0,N)-response );
    
    
end
toc;
%% Aggregate across participants
allw = cat(1,sub.w)';
allOptw = cat(1,sub.wOpt)';

wMean = mean(allw,2);
wSEM = std(allw,[],2)/sqrt(numel(subInd));

sMeanOpt = mean(cat(2,sub.wOpt),2);
%sMeanHeur = mean(cat(2,sub.slope_heur),2);

wMod_all = cat(3,sub.wMod);
sMeanMod = mean(wMod_all,3);


%% testing
for m=1:n_models
    [~, pval(m)] = ttest( allw(2,:),permute( wMod_all(m,2,:),[2 3 1]));
end

%% plot

figure(1);
width = 8;
height = 6;
LW = 1.2;
FS = 11;
clf;

hold on

barh(sMeanMod(:,2));
barh(n_models+1,wMean(2),'facecolor','g');
errorbar(wMean(2),n_models+1,wSEM(2),'horizontal','Color','k','LineWidth',LW);
set(gca, 'ytick', 1:n_models+1, 'yticklabel', [model_labels {'participants'}]);
xlabel('power dependence');


set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');

% figure name
figname = 'basic_N_regression_exponent';
set(gcf,'name',figname);


%% Print
filename =  fullfile('.\..\..\..\plots\exp1\supp_info\',[figname '.png']);
print(gcf, '-dpng', '-r400', filename);


