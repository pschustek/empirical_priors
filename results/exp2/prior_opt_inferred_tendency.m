%% Influence of the prior belief on responses for the optimal model and behavior 
% Fig. 5b
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\src\');

% Subselect subjects
subInd = [1:20 22:24];

% Load data
load('.\..\..\data\exp2_data.mat')

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

    pbConf = trials.prevConfBlockHeads;
    response = trials.confHeads;
    optRes = trials.optConfHeads;

    % Align with real block bias (symmetric pattern)
    block = 2*(trials.blockBias-0.5);
    response = block.*(response-0.5) + 0.5;
    optRes = block.*(optRes-0.5) + 0.5;
    mEv = block.*(mEv-0.5) + 0.5;
    pbConf = block.*(pbConf-0.5) + 0.5;

    %% Calculate bin values
    binIdx = discretize(pbConf,edges);     
    pbEvMask = binIdx==[1:nbins];
    
    for j=1:size(pbEvMask,2)
        sub(s).confDec(j) = mean(response(pbEvMask(:,j)),'omitnan');
        opt(s).confDec(j) = mean(optRes(pbEvMask(:,j)),'omitnan');
    end

    %% Fit bEv
    % very noisy
    lmod = fitlm(pbConf,optRes,'linear');
    opt(s).linConst = lmod.Coefficients.Estimate(1);
    opt(s).linSlope = lmod.Coefficients.Estimate(2);
    
    lmod = fitlm(pbConf,response,'linear');
    sub(s).linConst = lmod.Coefficients.Estimate(1);
    sub(s).linSlope = lmod.Coefficients.Estimate(2);

end

%% Aggregate across participants
allDecConf = cat(1,sub.confDec);
confMean = mean(allDecConf);
confSEM = std(allDecConf)/sqrt(numel(subInd));

confMeanOpt = mean(cat(1,opt.confDec));
confSEOpt = mean(cat(1,opt.confDec));

% Participants track close correlate of normative prior
X1 = vertcat(sub.confDec);
X2 = vertcat(opt.confDec);
[r, p] = corr(X1(:), X2(:), 'type', 'pearson');
fprintf('- [result] correlation with optimal: r = %.3f, p-value = %.3e\n', r, p);

% Less sensitivity to objective prior information
% Is the slope of the participants smaller than the optimal slope?
[p,~,stats] = signrank(vertcat(sub.linSlope),vertcat(opt.linSlope),'tail','left');
fprintf('- [result] smaller slopes, p-value = %.4f\n', p);

%% Plot
figure(1);
width = 8;
height = 6;
LW = 1.2;
FS = 11;
clf;

hold on
clear h

% q = linspace(0.2,1,100);
% plot(q,q, 'Color', [1 1 1]*0.7, 'LineWidth', 0.7, 'LineStyle', '-.');

h(2) = plot(binCenter,confMeanOpt,'Color', hsv2rgb([0 .7 0.9]) , 'LineWidth', LW+2);
h(1) = errorbar(binCenter,confMean,confSEM,'k', 'LineWidth', LW);

%xlim([0.5 blockLength+0.5]);
ylim([0.45 .85]);
xlabel('aligned inferred tendency');
ylabel('aligned confidence');

lh = legend(h,'exp','opt');
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
print(gcf, '-dpng', '-r400', '.\..\..\plots\exp2\prior_opt_inferred_tendency.png');