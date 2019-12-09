% Supp Fig. 8b
clear variables
close all

nu1 = 12;
nu2 = 9; % for heuristics model

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = [1:24]; 

% Load data
load('.\..\..\..\data\exp2_data.mat');

% Make fixed mEv bins
nbins = 4;
edges = [0 0.4545 0.6250 0.7500 1.0000];
binCenter = [diff(edges)/2 + edges(1:end-1)]';

%% Compute pattern
for s=subInd
    
    clear trials
    trials = trialData{s};
    
    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    
    % compute momentary evidence
    bEv = opt_inf.bEv( mEv.*N, N, nu1, nu2 );
    
    % Block evidence prior to present trial
    pbEv = [nan; bEv(1:end-1)];
    pbEv(1:trials.blockLength:numel(bEv)) = nan;
    
    response = trials.confHeads;
    optRes = trials.optConfHeads;
    
      
    % Align with real block bias
    block = 2*(trials.blockBias-0.5);
    response = block.*(response-0.5) + 0.5;
    optRes = block.*(optRes-0.5) + 0.5;
    mEv = block.*(mEv-0.5) + 0.5;
    pbEv = block.*(pbEv-0.5) + 0.5;
    
    % Remove first in-block trial
    idx = 1:trials.blockLength:numel(mEv);
    response(idx) = [];
    optRes(idx) = [];
    mEv(idx) = [];
    pbEv(idx) = [];
    
    %% Make groupings according to previous sample size
    % Mask for mEv
    binIdx = discretize(pbEv,edges);
    
    res = splitapply(@mean,response,binIdx);
    sub(s).confBin = res;
    
    res = splitapply(@mean,optRes,binIdx);
    opt(s).confBin = res;
        
    
    
end

%% Aggregate across participants
confMean = mean(cat(3,sub.confBin),3);
confSEM = std(cat(3,sub.confBin),[],3)/sqrt(numel(subInd));

confMeanOpt = mean(cat(3,opt.confBin),3);


%% Plot
figure(1);
clf;
width = 5.2;
height = width/1.333;
LW = 1.2;
FS = 7;
clf;
hold on
clear h

line([0 100],0.5*[1 1],'LineStyle','--','LineWidth',0.8,'Color',[1 1 1]*0.8);
line(50*[1 1],[0 1],'LineStyle','--','LineWidth',0.8,'Color',[1 1 1]*0.8);

% Model
blockLineStyle = {'-','-','-'};
colorMap = [hsv2rgb([.08 .5 0.85]); hsv2rgb([.33 .5 .85]) ; 0.59 0.82 0.59];

h(1) = plot(binCenter*100,confMeanOpt,'Color','b','LineStyle',blockLineStyle{1}, 'LineWidth', LW+2);

% Participants
blockLineStyle = {'-','-','-'};
colorMap = [hsv2rgb([.08 1 0.85]); hsv2rgb([.33 1 .85]); nan(1,3)];
h(2) =    errorbar(binCenter*100,confMean,confSEM,'Color','k','LineStyle',blockLineStyle{1}, 'LineWidth', LW);

xlim([0.2 .9]*100);
ylim([0.55 0.82]);

clear xL
xL = xlabel('aligned previous message (\%)  ','interpreter','latex','FontSize',11);
% set(xL,'Position',[53 xL.Position(2:end)]);
ylabel('aligned confidence');

legend(h, 'optimal', 'participants','location','northwest');

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');

%% save figure
figname = '.\..\..\..\plots\exp2\supp_info\prior_modulation_N_previous_message.png';
print(gcf, '-dpng', '-r400', figname);




