% Fig. 6a
clear variables
close all

%% Setup

nu1 = 12;
nu2 = 9; % for heuristics model

% Add path to auxiliary functions
addpath('.\..\..\src\');

% Subselect subjects
subInd = 1:24;

% Load data
load('.\..\..\data\exp2_data.mat');

%% Compute pattern
for s=subInd
    
    clear trials
    trials = trialData{s};
    
    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    
    response = trials.confHeads;
    optRes = trials.optConfHeads;
    
    %% Fit N
    setN = unique(N);
    x0 = 3;
    lb = 0.1;
    ub = 30;
    
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    
    % Normalize mEv
    nH = mEv .* N;
    mEv = 2*(mEv-0.5);
    
    ff = @(w,X) genlogf(X,w);
    
    w = nan(length(setN),2);
    R2exp = nan(length(setN),1);
    w_opt = nan(length(setN),2);
    w_avgN = nan(length(setN),2);
    R2opt = nan(length(setN),1);
    for j=1:length(setN)
        mask = setN(j)==N;
        
        
        % Participant
        fObj = @(w) -f_obj( ff(w,mEv(mask)), response(mask) );
        [w(j,:),~,exitflag] = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
        R2exp(j) = rsquared( ff(w(j,:), mEv(mask)), response(mask) );
        
        %         clf;
        %         hold on
        %         scatter(mEv(mask), response(mask));
        %         scatter(mEv(mask), ff(w(j,:), mEv(mask)),'r');
        
        % Optimal model
        fObj = @(w) -f_obj( ff(w,mEv(mask)), optRes(mask) );
        [w_opt(j,:),~,exitflag] = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
        R2opt(j) = rsquared( ff(w_opt(j,:), mEv(mask)), optRes(mask) );
        
        % Optimal model
        heurRes = opt_inf.all_approx( nH, N, trials.blockLength(1), nu1, nu2 );
        fObj = @(w) -f_obj( ff(w,mEv(mask)), heurRes(mask) );
        [w_heur(j,:),~,exitflag] = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
        R2heur(j) = rsquared( ff(w_opt(j,:), mEv(mask)), heurRes(mask) );
    end
    sub(s).w = w;
    sub(s).w_opt = w_opt;
    sub(s).w_heur = w_heur;
    sub(s).R2exp = R2exp;
    sub(s).R2opt = R2opt;
    sub(s).R2heur = R2opt;
    
end

%% Aggregate across participants
tmp = cat(1,sub.w);
allw = reshape(tmp(:,1), numel(setN),[]);
wMean = mean(allw,2);
wSEM = std(allw,[],2)/sqrt(numel(subInd));

tmp = cat(1,sub.w_opt);
allOptw = reshape(tmp(:,1), numel(setN),[]);
wOptMean = mean(allOptw,2);

tmp = cat(1,sub.w_heur);
allHeurw = reshape(tmp(:,1), numel(setN),[]);
wHeurMean = mean(allHeurw,2);

%% Positive slope on group level
tmp = repmat(setN,1,size(allw,2));

[rho, p] = corr(tmp(:),allw(:), 'type', 'pearson');
fprintf('- [result] correlation with sample size r = %.3f, p-value = %.3e\n', rho, p);

%% Plot
figure(1);
width = 5.2;
height = width/1.333;
LW = 1.2;
FS = 7;
clf;

hold on
clear h

colGrad = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*.5, ones(numel(q),1)*0.85]);
cmap = colGrad(0:0.02:1);
colormap(cmap);
patch([setN' nan],[wOptMean' nan],[linspace(0,1,numel(setN)) nan],'EdgeColor','interp','Marker','none','MarkerFaceColor','flat','LineWidth', LW+2);
%patch([setN' nan],[wHeurMean' nan],[linspace(0,1,numel(setN)) nan],'linestyle','--','EdgeColor','interp','Marker','none','MarkerFaceColor','flat','LineWidth', LW+2);


colGrad = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*1, ones(numel(q),1)*0.85]);
cmap = colGrad(0:0.02:1);
colormap(cmap);
h(1) = errorbar(setN,wMean,wSEM, 'Color', 'k', 'LineWidth', LW,'LineStyle','none');
patch([setN' nan],[wMean' nan],[linspace(0,1,numel(setN)) nan],'EdgeColor','interp','Marker','o',...
    'MarkerEdgeColor', 'k', 'LineWidth', .9,'MarkerFaceColor','flat','LineStyle', 'none', ...
    'MarkerSize', 4);

xlim([2 12]);
ylim([1 5]);
xlabel('sample size', 'FontSize', FS);
ylabel('slope confidence curve', 'FontSize', FS);

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');

%% Print
print(gcf, '-dpng', '-r400', '.\..\..\plots\exp2\prior_N_regression.png');