% Fig. 7b
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
    optRes = trials.optConfHeads;
    heurRes = opt_inf.all_approx( mEv.*N, N, trials.blockLength(1), nu1, nu2 );
    
    
    
    % Align with real block bias
    block = 2*(trials.blockBias-0.5);
    response = block.*(response-0.5) + 0.5;
    optRes = block.*(optRes-0.5) + 0.5;
    heurRes = block.*(heurRes-0.5) + 0.5;
    
    
    
    
    %% Median conditioned on in/block trial
    bMask = inBlock==1:blockLength;
    R = nan(numTrials/blockLength,blockLength);
    M = nan(numTrials/blockLength,blockLength);
    H = nan(numTrials/blockLength,blockLength);
    for j=1:blockLength
        R(:,j) = response(bMask(:,j));
        M(:,j) = optRes(bMask(:,j));
        H(:,j) = heurRes(bMask(:,j));
    end
    sub(s).res = R;
    sub(s).resOpt = M;
    sub(s).resHeur = H;
    sub(s).avg = mean(R);
    sub(s).avgOpt = mean(M);
    sub(s).avgHeur = mean(H);
    
end

%% Aggregate across participants
allMed = cat(1,sub.avg);
allOptw = cat(1,sub.avgOpt);
allHeurw = cat(1,sub.avgHeur);


wMean = mean(allMed,1);
wSEM = std(allMed,[],1)/sqrt(numel(subInd));

wMeanOpt = mean(allOptw);
wHeurOpt = mean(allHeurw);



%% Prior belief accumulation
% Four preditors
tmpY = allMed;
tmpX = repmat([1 2 3 4 5],numel(subInd),1);
% weight over trial index
mdl = fitlm(tmpX(:),tmpY(:),'linear');
w_pval = mdl.Coefficients.pValue(2);
% fprintf('>> Reported result in manuscript: p-value = %d\n', w_pval);
fprintf('- [result] p-value = %.3e\n', w_pval);

% Other test: Collapse bins (1,2)>1 and (4,5)>2
first = mean(allMed(:,[1 2]),2);
last = mean(allMed(:,[4 5]),2);
all_pval = signrank(last,first,'tail','right');

%% Plot
LW = 1.2;
FS = 11;
figure(1);
figname = 'prior_accum_block_length';


width = 8;
height = 6.5;
clf;

legLabel = {'exp','opt'};

hold on
clear h

h(1) = errorbar(1:blockLength, wMean, wSEM, 'Color','k', 'LineWidth', LW);
h(2) = plot(1:blockLength,wMeanOpt, 'Color', hsv2rgb([0 .7 0.9]) , 'LineWidth', LW+2);

xlim([0.5 blockLength+0.5]);
ylim([0.6 0.75]);

xlabel('trial in block');
ylabel('aligned confidence');

lh = legend(h,legLabel,'Location','northwest');
set(lh,'Box','off','Interpreter','latex');

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',1:blockLength);

%% Print
filename = fullfile( '.\..\..\plots\exp2', [figname '.png']);
print(gcf, '-dpng', '-r400', filename);