%% Calibration curve
% Fig. S1b
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = [1:20 22:24];

% Load data
load('.\..\..\..\data\exp2_data.mat');

% Prepare data
nbins = 8;
% Pool all responses together
accum = [];
for s=subInd
    accum = [accum; trialData{s}.confDec];
end
% % With respect to equally filled bins
qtls = [0:nbins]/nbins;
edges = quantile(accum,qtls(2:end-1));
edges = [0.5 edges 1];
binCenter = [diff(edges)/2 + edges(1:end-1)]';

%% Compute
for s=subInd

    clear trials
    trials = trialData{s};
    numTrials = size(trials.num,1);
    
    C = trials.confDec;
    
    % Objective probability of being correct for decision [1,0]
    objDecCorr = [trials.optConfHeads 1-trials.optConfHeads];
    
    % Probability correct for participant's response
    probDecCorr = sum(objDecCorr.*[trials.decision ~trials.decision],2);
    
    optC = abs(trials.optConfHeads-0.5) + 0.5;
    optCorr = trials.optDecHeads == trials.coinBias;
    
    % Participants
    binIdx = discretize(C,edges);     
    mask = binIdx==[1:nbins];
    for j=1:size(mask,2)
        sub(s).probCorr(j) = mean(probDecCorr(mask(:,j)),'omitnan');
    end
    
    % Optimal model
    binIdx = discretize(optC,edges);     
    mask = binIdx==[1:nbins];
    for j=1:size(mask,2)
        opt(s).probCorr(j) = mean(optCorr(mask(:,j)),'omitnan');
    end
    
end

%% Aggregate across participants
Q = cat(1,sub.probCorr);
Qopt = cat(1,opt.probCorr);

allMean = mean(Q,'omitnan');
allSEM = std(Q,'omitnan')/sqrt(size(Q,1));

% Tests of deviation from unity line
for j=1:size(Q,2)
    tmp = Q(:,j);
    pval.above(j) = signrank(tmp(~isnan(tmp)),binCenter(j),'tail','right');  
    pval.below(j) = signrank(tmp(~isnan(tmp)),binCenter(j),'tail','left');
end

% Group level calibratedness: Correlation with model
mask = ~isnan(Q);
[rho,pcal] = corr(Q(mask),Qopt(mask));
fprintf('- [result] correlation with model, r = %.3f, p = %.3e\n', rho, pcal);

%% Plot
width = 7.5;
height = 6;
LW = .9;
FS = 9;
figure(1)
clf;

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

hold on

line([0.5 1],[0.5 1],'Color',[1 1 1]*0.5,'LineWidth',1,'LineStyle','--');
errorbar(binCenter,allMean,allSEM,'Color','k','LineWidth',LW,'LineStyle','-','CapSize',3);

for j=1:6
    if pval.above(j)<=0.05 & pval.above(j)>0.01
        text(binCenter(j),allMean(j)+0.07,'$\ast$','Interpreter','latex','FontSize', FS-1,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    if pval.above(j)<=0.01 
        text(binCenter(j),allMean(j)+0.07,'$\ast\ast$','Interpreter','latex','FontSize', FS-1,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

for j=7:length(binCenter)
    if pval.below(j)<=0.05 & pval.below(j)>0.01
        text(binCenter(j),allMean(j)-0.07,'$\ast$','Interpreter','latex','FontSize', FS-1,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    if pval.below(j)<=0.01 
        text(binCenter(j),allMean(j)-0.07,'$\ast\ast$','Interpreter','latex','FontSize', FS-1,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

xlim([0.5 1]);
ylim([0.5 1]);

xlabel('confidence in decision', 'FontSize', FS, 'FontName', 'Times');
ylabel('probability correct', 'FontSize', FS, 'FontName', 'Times');

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','YTick',0.5:0.1:1,'XTick',0.5:0.1:1);


%% Print
print(gcf, '-dpng', '-r400', '.\..\..\..\plots\exp2\supp_info\si_prior_confidence_calibration.png');