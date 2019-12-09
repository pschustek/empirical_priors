% Supp Fig. 9
clear variables
close all

nu1 = 14;
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
    blockLength = unique(trials.blockLength);
    numTrials = size(trials.num,1);
    
    response = trials.confHeads;
    optRes = trials.optConfHeads;
    inBlock = trials.trialInBlock;
    
    
    % Symmetrize
    E = 2*(mEv-0.5);
    
    % Index corresponds to previous trial
    Xdat = zeros(numTrials,blockLength-1);
    numP = inBlock-1;
    for j=1:4
        Xdat(:,j) = [nan(j,1);E(1:end-j)]; % shift by j trials
        % x = E(any(inBlock==[1:j],2));
        % Xdat(any(inBlock==[5-(j-1):5],2),j) = x;
    end
    %  Xdat(inBlock==1,:) = [];
    %  numP(inBlock==1,:) = [];
    
    prev(s).wExp = fit_prior_weights_fixed( fliplr(Xdat), response, numP );
    prev(s).wMod = fit_prior_weights_fixed( fliplr(Xdat), optRes, numP );
end

%% Aggregate across participants
wExp = cat(3,prev.wExp);
wMod = cat(3,prev.wMod);

modMean = nanmean(wMod,3);
expMean = nanmean(wExp,3);
expSE = nanstd(wExp,[],3)/sqrt(numel(subInd));

for i=1:5
    for j=1:4
        [~, pval(i,j)] = ttest(squeeze(wExp(i,j,:)));
        if j>=i
            fprintf('p-value for trial lag t-%d at block position %d: %f\n', j, i, pval(i,j));
        end
    end
end

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
width = 2*8;
height = 2*6.5;
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
%subplot('Position',[xedge(1) yedge(3) xe(1) ye(3)]);
for i=1:5
    % weights corresponding to previous block trials
    subplot(5,2,2*i-1); title(['trial position ' num2str(i)]);
    hold on
    if i<5
        idx = 1:5-i;
        h(1) = plot([i+1:5]', modMean(i,idx), 'Color', colMod(1,:), 'LineWidth', LW+2, 'Marker', '.');
        h(2) = errorbar([i+1:5]',expMean(i,idx),expSE(i,idx), 'Color', colExp(1,:), 'LineWidth', LW);
    end
    axis tight; xlim([1.4 5.6]);
    plot(xlim,[0 0], 'color','k', 'linewidth',0.5);
    
    
    
    if i==3
        ylabel('weight of sample proportion ', 'Interpreter','latex');
    end
    if i==5
        xlabel('previous block trial position');
    else
        set(gca, 'xticklabel','');
    end
    
    % weights corresponding to current block trials
    subplot(5,2,2*i);
    hold on
    if  i>1
        idx = 6-i:4;
        h(1) = plot([1:i-1]', modMean(i,idx), 'Color', colMod(1,:), 'LineWidth', LW+2, 'Marker', '.');
        h(2) = errorbar([1:i-1]',expMean(i,idx),expSE(i,idx), 'Color', colExp(1,:), 'LineWidth', LW);
    end
    axis tight; xlim([0.4 4.6]);
    plot(xlim,[0 0], 'color','k', 'linewidth',0.5);
    
    if i==1
        lh = legend(h,'opt','exp','Location','northeast');
        set(lh,'Box','off','Interpreter','latex');
    end
    if i==5
        xlabel('current block trial position');
    else
        set(gca, 'xticklabel','');
    end
    set(gca, 'yticklabel','');
    
end
sameaxis;


% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',1:4);

%% Print
print(gcf, '-dpng', '-r400', '.\..\..\plots\exp2\prior_weight_previous_trials_spillover.png');