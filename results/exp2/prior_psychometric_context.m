%% Influence of the real block tendency on behavior
% Fig. 5a
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\src\');

% Subselect subjects
subInd = [1:20 22:24];

% Load data
load('.\..\..\data\exp2_data.mat');

% Prepare data
nbins = 9;
% Pool all responses together
accum = [];
for s=subInd
    accum = [accum; trialData{s}.meanEvidence];
end
% % With respect to equally filled bins
qtls = [0:nbins]/nbins;
edges = quantile(accum,qtls(2:end-1));
edges = [0 edges 1];
binCenter = [diff(edges)/2 + edges(1:end-1)]';

%% Prepare data
for s=subInd

    clear trials    
    trials = trialData{s};

    mEv = trials.meanEvidence;
    response = trials.confHeads;
    optRes = trials.optConfHeads;

    binIdx = discretize(mEv,edges);     
    [G,~,bb] = findgroups(binIdx, trials.blockBias);
    res = splitapply(@mean,response,G);
    sub(s).confBin = [res(bb==1) res(bb==0)];
    
    res = splitapply(@mean,optRes,G);
    opt(s).confBin = [res(bb==1) res(bb==0)];

end

%% Aggregate across participants
confMean = mean(cat(3,sub.confBin),3);
confSEM = std(cat(3,sub.confBin),[],3)/sqrt(numel(subInd));

confMeanOpt = mean(cat(3,opt.confBin),3);

%% Plot
figure(1);
width = 8;
height = 6;
LW = 1.2;
FS = 11;
clf;

hold on

h(1) = plot(binCenter*100,confMeanOpt(:,1), 'Color', hsv2rgb([.5 .5 0.9]), 'LineWidth', LW+2);
h(2) = plot(binCenter*100,confMeanOpt(:,2), 'Color', hsv2rgb([.86 .5 0.9]), 'LineWidth', LW+2);

q = linspace(0,100,200);
plot(q,q/100, 'Color', [1 1 1]*0.8, 'LineWidth', 0.8, 'LineStyle', '--');

errorbar(binCenter*100,confMean(:,1),confSEM(:,1),'Color', hsv2rgb([.5 1 0.9]), 'LineWidth', LW);
errorbar(binCenter*100,confMean(:,2),confSEM(:,2),'Color', hsv2rgb([.86 1 0.9]), 'LineWidth', LW);

%xlim([min(edges)-0.05 max(edges)+0.05]);
%ylim([-0.2 1]);
xlabel('proportion of blue samples (\%)','Interpreter','latex');
ylabel('confidence blue majority','Interpreter','latex');

% Custom legend
xl = xlim;
yl = ylim;
xpos = 7;
ypos = 0.8;
rectangle('Position',[xpos ypos diff(xl)/10 diff(yl)/40],'EdgeColor','none','FaceColor',hsv2rgb([.5 .5 0.9]));
ypos = 0.7;
rectangle('Position',[xpos ypos diff(xl)/10 diff(yl)/40],'EdgeColor','none','FaceColor',hsv2rgb([.86 1 0.9]));

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
        'XMinorTick', 'off', 'YMinorTick', 'off', 'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');

%% Print
print(gcf, '-dpng', '-r400', '.\..\..\plots\exp2\prior_psychometric_context.png');