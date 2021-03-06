% Fig. 7a
clear variables
close all

nu1 = 12;
nu2 = 9; % for heuristics model

%% Setup
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
    blockLength = unique(trials.blockLength);
    numTrials = size(trials.num,1);
    
    response = trials.confHeads;
    optRes = trials.optConfHeads;
    inBlock = trials.trialInBlock;
    
    % Dependent variable
    mask = any(inBlock==[2:5],2);
    expY = response(mask);
    modelY = optRes(mask);
    
    % Symmetrize
    E = 2*(mEv-0.5);
    
    
    % Index corresponds to previous trial
    Xdat = zeros(numTrials,blockLength-1);
    numP = inBlock-1;
    for j=1:4
        x = E(any(inBlock==[1:j],2));
        Xdat(any(inBlock==[5-(j-1):5],2),j) = x;
    end
    Xdat(inBlock==1,:) = [];
    numP(inBlock==1,:) = [];
    
    prev(s).wExp = fit_prior_weights( fliplr(Xdat), expY, numP );
    prev(s).wOpt = fit_prior_weights( fliplr(Xdat), modelY, numP );
end

%% Aggregate across participants
wExp = cat(3,prev.wExp);
wOpt = cat(3,prev.wOpt);

optMean = nanmean(wOpt,3);

expMean = nanmean(wExp,3);
expSE = nanstd(wExp,[],3)/sqrt(numel(subInd));

%% Regression of weight over trials index
% Two preditors
tmpY = squeeze(wExp(2,1:2,:));
tmpX = repmat([1 2]',1,numel(subInd));
% weight over trial index
mdl = fitlm(tmpX(:),tmpY(:),'linear');
w_pval(1) = mdl.Coefficients.pValue(2);
tmpYall{1} = tmpY;

% Three preditors
tmpY = squeeze(wExp(3,1:3,:));
tmpX = repmat([1 2 3]',1,numel(subInd));
% weight over trial index
mdl = fitlm(tmpX(:),tmpY(:),'linear');
w_pval(2) = mdl.Coefficients.pValue(2);
tmpYall{2} = tmpY;


% Four preditors
tmpY = squeeze(wExp(4,1:4,:));
tmpX = repmat([1 2 3 4]',1,numel(subInd));
% weight over trial index
mdl = fitlm(tmpX(:),tmpY(:),'linear');
w_pval(3) = mdl.Coefficients.pValue(2);
tmpYall{3} = tmpY;

fprintf('- [result] p-values of regressions for 2-4 previous trials = %.3f, %.3f, %.3f\n', w_pval);
for j=1:3
    [~, pval_last_two(j)] = ttest(tmpYall{j}(end-1,:), tmpYall{j}(end,:)); % comparing weight for trials t-2 and t-1
    fprintf('- [result] p-values for comparing weights at t-2 and t-1 for trial position %d: = %.3f\n', j+2, pval_last_two(j));
    tmpY_last_two{j} = tmpYall{j}(end-1:end,:);
    
end
tmpY_altogether = [tmpY_last_two{:}];
[~, pval_last_two_all] = ttest(tmpY_altogether(1,:), tmpY_altogether(2,:)); % comparing weight for trials t-2 and t-1
fprintf('- [result] p-values for comparing weights at t-2 and t-1 altogether: = %.3f\n', pval_last_two_all);


%% Plot
colMod = repmat(hsv2rgb([0 .7 0.9]),3,1);
colExp = repmat([0 0 0],3,1);
colFit = repmat([0.5 0.5 1],3,1);
ybnd = [0.5 1.5];

nXpanels = 1;
nYpanels = 3;
x_lmargin = 0.17;
x_rmargin = 0.05;
x_relWidth = [1];
x_spacing = 0.0;
y_lmargin = 0.2;
y_umargin = 0.04;
y_relWidth = [0.333 0.333 0.333];
y_spacing = 0.08;

% x and y extent
xe = (1-(x_lmargin+x_rmargin+2*x_spacing))*x_relWidth;
ye = (1-(y_lmargin+y_umargin+2*y_spacing))*y_relWidth;
% Left edge
sep = [x_lmargin];  % all separations
% assert(sum(sep)==1);
sep = cumsum(sep);
xedge = [sep(1)];
% Lower edge
sep = [y_lmargin ye(1) y_spacing ye(2) y_spacing ye(3) y_umargin];  % all separations
% assert(abs(sum(sep)-1)<=0.01);
sep = cumsum(sep);
yedge = [sep(1) sep(3) sep(5)];

figure(1);
figname = 'prior_weight_previous_trials';
width = 8;
height = 6.5;
LW = 1.2;
FS = 11;
clf;
clear h

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% First
subplot('Position',[xedge(1) yedge(3) xe(1) ye(3)]);
hold on
h(1) = plot([1:4]', optMean(2,:), 'Color', colMod(1,:), 'LineWidth', LW+2, 'Marker', '.');
h(2) = errorbar([1:4]',expMean(2,:),expSE(2,:), 'Color', colExp(1,:), 'LineWidth', LW);


xlim([0.5 4.5]);
ylim(ybnd);

lh = legend(h,'opt','exp','Location','best');
set(lh,'Box','off','Interpreter','latex');

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',1:4,'XTickLabel',{});

% Second
subplot('Position',[xedge(1) yedge(2) xe(1) ye(2)]);
hold on
h(1) = plot([1:4]', optMean(3,:), 'Color', colMod(2,:), 'LineWidth', LW+2, 'Marker', '.');
h(2) = errorbar([1:4]',expMean(3,:),expSE(3,:), 'Color', colExp(2,:), 'LineWidth', LW);
xlim([0.5 4.5]);
ylim(ybnd);

set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',1:4,'XTickLabel',{});

ylabel('weight of sample proportion ', 'Interpreter','latex');

% Third
subplot('Position',[xedge(1) yedge(1) xe(1) ye(1)]);
hold on
h(1) = plot([1:4]', optMean(4,:), 'Color', colMod(3,:), 'LineWidth', LW+2, 'Marker', '.');
h(2) = errorbar([1:4]',expMean(4,:),expSE(4,:), 'Color', colExp(3,:), 'LineWidth', LW);
xlim([0.5 4.5]);
ylim(ybnd);

xlabel('trial in block');

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',1:4);

%% Print
filename = fullfile( '.\..\..\plots\exp2', [figname '.png']);
print(gcf, '-dpng', '-r400', filename);