% Fig. S4d
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
    numTrials = size(trials.num,1);
    blockInd = ceil([1:numTrials]'/unique(trials.blockLength));

    % Block evidence prior to present trial
    pmEv = [nan; mEv(1:end-1)];
    pmEv(1:trials.blockLength:numel(mEv)) = nan;
    pN = [nan; trials.sampleSize(1:end-1)];
    pN(1:trials.blockLength:numel(mEv)) = nan;

    response = trials.confHeads;
    
    % Fitted model pHT_pHT_sig
    w = ml_prior(s).pHT_pHT_sig.weights;
    fitRes = pHT_pHT_sig(w,[mEv N blockInd]);

    % Align with real block bias
    block = 2*(trials.blockBias-0.5);
    response = block.*(response-0.5) + 0.5;
    mEv = block.*(mEv-0.5) + 0.5;
    pmEv = block.*(pmEv-0.5) + 0.5;
    fitRes = block.*(fitRes-0.5) + 0.5;
    
    % Remove first in-block trial
    idx = 1:trials.blockLength:numel(mEv);
    response(idx) = [];
    mEv(idx) = [];
    pmEv(idx) = [];
    pN(idx) = [];
    fitRes(idx) = [];
    
    % Mask for pmEv 
    binIdx = discretize(pmEv,edges);
    
    few = pN <= quantile(pN,0.5);
    many = pN > quantile(pN,0.5);
    groupMask = logical([few many]);    
    allConditions = any(groupMask,2);   
    
    [G,~,bb] = findgroups(binIdx(allConditions), groupMask(allConditions,1));
    res = splitapply(@mean,response(allConditions),G);
    sub(s).confBin = [res(bb==1) res(bb==0)];
    
    res = splitapply(@mean,fitRes(allConditions),G);
    fit(s).confBin = [res(bb==1) res(bb==0)];

end

%% Aggregate across participants
confMean = mean(cat(3,sub.confBin),3);
confSEM = std(cat(3,sub.confBin),[],3)/sqrt(numel(subInd));

confMeanFit = mean(cat(3,fit.confBin),3);

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

% Fitted model
blockLineStyle = {'-','-','-'};
colorMap = [hsv2rgb([.08 .5 0.85]); hsv2rgb([.33 .5 .85]); nan(1,3)];
legIdx = 3:4;
for b=1:size(groupMask,2)
    h(legIdx(b)) = plot(binCenter*100,confMeanFit(:,b),'Color',colorMap(b,:),'LineStyle',blockLineStyle{b}, 'LineWidth', LW+3);
end

% Participants
blockLineStyle = {'-','-','-'};
colorMap = [hsv2rgb([.08 1 0.85]); hsv2rgb([.33 1 .85]); nan(1,3)];
legIdx = 1:2;
for b=1:size(groupMask,2)
    h(legIdx(b)) = errorbar(binCenter*100,confMean(:,b),confSEM(:,b),'Color',colorMap(b,:),'LineStyle',blockLineStyle{b}, 'LineWidth', LW);
end

xlim([0.2 .9]*100);
ylim([0.55 0.76]);

xlabel('aligned proportion previous sample (\%)  ','interpreter','latex');
ylabel('aligned confidence');

% Custom legend
xl = xlim;
yl = ylim;
xpos = 24;
ypos = 0.715;
rectangle('Position',[xpos ypos diff(xl)/8 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.33 1 .85]));
ypos = 0.69;
rectangle('Position',[xpos ypos diff(xl)/8 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.08 1 0.85]));

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],... 
    'XMinorTick', 'off', 'YMinorTick', 'off', 'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');


%% Print
print(gcf, '-dpng', '-r400', '.\..\..\..\plots\exp2\supp_info\si_prior_modulation_N_previous_trial.png');