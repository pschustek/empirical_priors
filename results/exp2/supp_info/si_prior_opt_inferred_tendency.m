%% Influence of the prior belief on responses for the optimal model and behavior 
% Fig. S4b
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = [1:20 22:24];

% Load data
load('.\..\..\..\data\exp2_data.mat');
load('.\..\..\..\models\exp2\prior_ml.mat');

% Prepare data
nbins = 7;
% Pool all responses together
accum = [];
for s=subInd
    block = 2*(trialData{s}.blockBias-0.5);
    accum = [accum; block.*(trialData{s}.prevConfBlockHeads-0.5)+0.5];
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
    N = trials.sampleSize;
    numTrials = size(trials.num,1);
    blockInd = ceil([1:numTrials]'/unique(trials.blockLength));

    pbConf = trials.prevConfBlockHeads;
    response = trials.confHeads;
    
    % Fitted model pHT_pHT_sig
    w = ml_prior(s).pHT_pHT_sig.weights;
    fitRes = pHT_pHT_sig(w,[mEv N blockInd]);

    % Align with real block bias (symmetric pattern)
    block = 2*(trials.blockBias-0.5);
    response = block.*(response-0.5) + 0.5;
    mEv = block.*(mEv-0.5) + 0.5;
    pbConf = block.*(pbConf-0.5) + 0.5;
    fitRes = block.*(fitRes-0.5) + 0.5;

    %% Calculate bin values
    binIdx = discretize(pbConf,edges);     
    pbEvMask = binIdx==[1:nbins];
    
    for j=1:size(pbEvMask,2)
        sub(s).confDec(j) = mean(response(pbEvMask(:,j)),'omitnan');
        fit(s).confDec(j) = mean(fitRes(pbEvMask(:,j)),'omitnan');
    end
end

%% Aggregate across participants
allDecConf = cat(1,sub.confDec);
confMean = mean(allDecConf);
confSEM = std(allDecConf)/sqrt(numel(subInd));

confMeanFit = mean(cat(1,fit.confDec));

%% Plot
figure(1);
width = 5.2;
height = width/1.333;
LW = 1.2;
FS = 7;
clf;

hold on
clear h

h(2) = plot(binCenter,confMeanFit,'Color', hsv2rgb([0 .7 0.9]), 'LineWidth', LW+2, 'LineStyle', '-');
h(1) = errorbar(binCenter,confMean,confSEM,'k', 'LineWidth', LW);

%xlim([0.5 blockLength+0.5]);
ylim([0.48 0.8]);
xlabel('aligned inferred tendency');
ylabel('aligned confidence');

lh = legend(h,'exp','fit');
set(lh,'Interpreter','latex','Box','off','Location','best');

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
print(gcf, '-dpng', '-r400', '.\..\..\..\plots\exp2\supp_info\si_prior_opt_inferred_tendency.png');