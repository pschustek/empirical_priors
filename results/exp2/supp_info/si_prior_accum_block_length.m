% Fig. S4f
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

%% Compute patttern
for s=subInd

    clear trials    
    trials = trialData{s};

    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    blockLength = unique(trials.blockLength);
    numTrials = size(trials.num,1);
    inBlock = trials.trialInBlock;
    blockInd = ceil([1:numTrials]'/unique(trials.blockLength));

    response = trials.confHeads;

    % Fitted model pHT_pHT_sig
    w = ml_prior(s).pHT_pHT_sig.weights;
    fitRes = pHT_pHT_sig(w,[mEv N blockInd]);

    % Align with real block bias
    block = 2*(trials.blockBias-0.5);
    response = block.*(response-0.5) + 0.5;
    mEv = block.*(mEv-0.5) + 0.5;
    fitRes = block.*(fitRes-0.5) + 0.5;

    %% Median conditioned on in/block trial
    bMask = inBlock==1:blockLength;
    R = nan(numTrials/blockLength,blockLength);
    F = nan(numTrials/blockLength,blockLength);
    for j=1:blockLength
        R(:,j) = response(bMask(:,j));
        F(:,j) = fitRes(bMask(:,j));
    end
    
    sub(s).res = R;
    sub(s).avg = mean(R);
    sub(s).fit = mean(F);
end

%% Aggregate across participants
allMed = cat(1,sub.avg);
allFitw = cat(1,sub.fit);

wMean = mean(allMed,1);
wSEM = std(allMed,[],1)/sqrt(numel(subInd));

wMeanFit = mean(allFitw);

%% Plot
width = 5.2;
height = width/1.333;
LW = 1.2;
FS = 7;
clf;

hold on

plot(1:blockLength,wMeanFit, 'Color', hsv2rgb([0 .7 0.9]), 'LineWidth', LW+2);
errorbar(1:blockLength, wMean, wSEM, 'Color','k', 'LineWidth', LW);


xlim([0.5 blockLength+0.5]);
ylim([0.6 0.71]);

xlabel('trial in block');
ylabel('aligned confidence');


set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
        'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',1:blockLength, 'YTick', 0.6:0.05:0.7);

%% Print
print(gcf, '-dpng', '-r400', '.\..\..\..\plots\exp2\supp_info\si_prior_accum_block_length.png');