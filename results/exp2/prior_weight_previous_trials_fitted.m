% Fig. 8e
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\src\');

% Subselect subjects
subInd = 1:24;

% Load data
load('.\..\..\data\exp2_data.mat');
load('.\..\..\models\exp2\prior_cvll.mat');

%% Prepare data
for s=subInd

    clear trials    
    trials = trialData{s};

    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    blockLength = unique(trials.blockLength);
    numTrials = size(trials.num,1);

    response = trials.confHeads;
    inBlock = trials.trialInBlock;
    blockIdx = ceil([1:numTrials]/blockLength)';
    
    % Fitted model pHT_pHT_sig
   fitRes =   cvll_prior(s).pHT_pHT_sig.Y; % model

    % Dependent variable
    mask = any(inBlock==[2:5],2);
    expY = response(mask);
    fitY = fitRes(mask);
    
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
    prev(s).wFit = fit_prior_weights( fliplr(Xdat), fitY, numP );

end
    
%% Aggregate across participants
wExp = cat(3,prev.wExp);
wFit = cat(3,prev.wFit);

fitMean = nanmean(wFit,3);
expMean = nanmean(wExp,3);
expSE = nanstd(wExp,[],3)/sqrt(numel(subInd));
    
%% Plot
colMod = repmat(hsv2rgb([0 .7 0.9]),3,1);
colExp = repmat([0 0 0],3,1);
colFit = repmat(hsv2rgb([0 .7 0.9]),3,1);
ybnd = [0.5 1.5];

nXpanels = 1;
nYpanels = 3;
x_lmargin = 0.2;
x_rmargin = 0.08;
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
width = 5.2;
height = width/1.333;
LW = 1.2;
FS = 7;
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
h(2) = plot([1:4]', fitMean(2,:), 'Color', colFit(1,:), 'LineWidth', LW+2, 'Marker', '.');
h(1) = errorbar([1:4]',expMean(2,:),expSE(2,:), 'Color', colExp(1,:), 'LineWidth', LW, 'Marker', 'none');
xlim([0.5 4.5]);
ylim(ybnd);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',1:4,'XTickLabel',{});

% Second
subplot('Position',[xedge(1) yedge(2) xe(1) ye(2)]);
hold on
h(2) = plot([1:4]', fitMean(3,:), 'Color', colFit(2,:), 'LineWidth', LW+2, 'Marker', '.');
h(1) = errorbar([1:4]',expMean(3,:),expSE(3,:), 'Color', colExp(2,:), 'LineWidth', LW, 'Marker', 'none');
xlim([0.5 4.5]);
ylim(ybnd);

set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',1:4,'XTickLabel',{});

ylabel('weight of sample proportion ', 'Interpreter','latex');

% Third
subplot('Position',[xedge(1) yedge(1) xe(1) ye(1)]);
hold on
h(2) = plot([1:4]', fitMean(4,:), 'Color', colFit(3,:), 'LineWidth', LW+2, 'Marker', '.');
h(1) = errorbar([1:4]',expMean(4,:),expSE(4,:), 'Color', colExp(3,:), 'LineWidth', LW, 'Marker', 'none');
xlim([0.5 4.5]);
ylim(ybnd);

xlabel('trial in block');

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',1:4);

%% Print
print(gcf, '-dpng', '-r400', '.\..\..\plots\exp2\prior_weight_previous_trials_fitted.png');