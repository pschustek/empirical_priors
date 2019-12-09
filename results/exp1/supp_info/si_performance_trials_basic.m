%% Supp Fig 6a; Evolution of performance index over trials
clear variables
close all

%% Prepare
% Data
load('.\..\..\..\data\exp1_data.mat');
addpath ././../../src

%% Prepare data
subInd =1:24;
nbins = 10;
edges = linspace(1,280,nbins+1);
binCenter = [diff(edges)/2 + edges(1:end-1)]';

%%
nfirst = 50;
devFirst = nan(max(subInd),nfirst);
perfFirst = nan(max(subInd),nfirst);
for s=subInd
    clear trials
    trials = trialData{s};
    
    deviation = trials.optConfHeads - trials.confHeads;
    numTrials = size(trials.num,1);
    
    binIdx = discretize([1:numTrials]',edges);     
    mask = binIdx==[1:nbins];
    
    for j=1:size(mask,2)
        D = deviation(mask(:,j));
        dev(s).mean(j) = mean(D);
        dev(s).median(j) = median(D);
        dev(s).std(j) = std(D);
        dev(s).SE(j) = std(D)/sqrt(length(D));
        
        P = abs(D);
        perf(s).mean(j) = mean(P);
        perf(s).median(j) = median(P);
        perf(s).std(j) = std(P);
        perf(s).SE(j) = std(P)/sqrt(length(P));
    end
    
    devFirst(s,:) = deviation(1:nfirst)';
    perfFirst(s,:) = abs(deviation(1:nfirst)');
end

%% Aggregate across participants
Qdev = cat(1,dev.mean);
devMean = mean(Qdev);
devSEM = std(Qdev)/sqrt(size(Qdev,1));

Qperf = cat(1,perf.mean);
perfMean = mean(Qperf);
perfSEM = std(Qperf)/sqrt(size(Qperf,1));

% Test if different from average response across all participants
j = 0;
R = nan(numel(subInd),1);
for s=subInd
    j = j + 1;
    deviation = trialData{s}.optConfHeads-trialData{s}.confHeads;
    numTrials = size(trialData{s}.num,1);
    D1(j) = mean(deviation(1:floor(numTrials/2)));
    D2(j) = mean(deviation(floor(numTrials/2)+1:end));
    
    P = abs(deviation);
    P1(j) = mean(P(1:floor(numTrials/2)));
    P2(j) = mean(P(floor(numTrials/2)+1:end));
end
signrank(D1,D2);
signrank(P1,P2);


%% Plot performance
width = 16;
height = 12;
LW = 1.5;
FS = 16;
figure(2)
clf;

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

hold on

errorbar(binCenter,perfMean,perfSEM,'Color','k','LineWidth',LW,'LineStyle','-');

xlim([0 281]);
ylim([0.08 .15]);
xlabel('trial', 'FontSize', FS, 'FontName', 'Times');
ylabel('absolute deviation from optimal', 'FontSize', FS, 'FontName', 'Times');


% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');

% Inset
inset = axes('Position',[0.25 0.7 0.35 0.2]);
hold on
SEM = std(perfFirst,'omitnan')/sqrt(numel(subInd));
plot(1:nfirst,mean(perfFirst,'omitnan'),'Color','k','LineWidth',LW,'LineStyle','-');
plot(1:nfirst,mean(perfFirst,'omitnan')+SEM,'Color','r','LineWidth',1,'LineStyle','-.');
plot(1:nfirst,mean(perfFirst,'omitnan')-SEM,'Color','r','LineWidth',1,'LineStyle','-.');

xlim([1 nfirst]);
ylim([0.05 .15]);

set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times','XTick',[1 25 50]);

% Print
print(gcf, '-dpng', '-r400', '.\..\..\..\plots\exp1\supp_info\performance_trials_basic.png');
