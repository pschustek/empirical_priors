%% Supp Fig 6b: Evolution of performance index over trials across sessions
clear variables
close all

%% Prepare
subInd = 1:24;

% Basic task data
%load('.\..\..\..\data\exp1_data.mat');
addpath('.\..\..\..\src\');

% Prior learning task data
load('.\..\..\..\data\exp2_data.mat');
for s=subInd
    trials = trialData{s};
    idx = [trials.num; 1];
    maxTrialsSession(s,1:2) = trials.num(diff(idx)~=1);
    devSub{s,1} = trials.confHeads(1:maxTrialsSession(s,1)) - trials.optConfHeads(1:maxTrialsSession(s,1));
    devSub{s,2} = trials.confHeads(maxTrialsSession(s,1)+1:end) - trials.optConfHeads(maxTrialsSession(s,1)+1:end);
    rtSub{s,1} = trials.responseTime(1:maxTrialsSession(s,1));
    rtSub{s,2} = trials.responseTime(maxTrialsSession(s,1)+1:end);
end

% Number of bin for each of the three sessions
nbins = [4 4];

%%
nfirst = 50;
devFirst = nan(max(subInd),nfirst);
perfFirst = nan(max(subInd),nfirst);
rtFirst = nan(max(subInd),nfirst);
for s=subInd
    
    % Index runs over all bins of all sessions together
    k = 0;
    binCenter = [];

    % Experimental sessions
    for b=1:2
        edges = linspace(1,maxTrialsSession(s,b),nbins(b)+1);
        binCenter = [binCenter; [diff(edges)/2 + edges(1:end-1)]'];
        
        binIdx = discretize([1:maxTrialsSession(s,b)]',edges);
        mask = binIdx==[1:nbins(b)];
        
        % Bins (within session)
        for j=1:size(mask,2)
            k = k + 1;
            D = devSub{s,b}(mask(:,j));
            dev(s).mean(k) = mean(D);
            dev(s).median(k) = median(D);
            dev(s).std(k) = std(D);
            dev(s).SE(k) = std(D)/sqrt(length(D));
            
            P = abs(D);
            perf(s).mean(k) = mean(P);
            perf(s).median(k) = median(P);
            perf(s).std(k) = std(P);
            perf(s).SE(k) = std(P)/sqrt(length(P));
            
            R = rtSub{s,b}(mask(:,j));
            RT(s).mean(k) = mean(R);
            RT(s).median(k) = median(R);
            RT(s).std(k) = std(R);
            RT(s).SE(k) = std(R)/sqrt(length(R));
        end
    end
    
    % First trials of first prior learning session (b=2)
    devFirst(s,:) = devSub{s,2}(1:nfirst)';
    perfFirst(s,:) = abs(devSub{s,2}(1:nfirst)');
    rtFirst(s,:) = rtSub{s,2}(1:nfirst)';
end

%% Aggregate across participants
Qdev = cat(1,dev.mean);
devMean = mean(Qdev);
devSEM = std(Qdev)/sqrt(size(Qdev,1));

Qperf = cat(1,perf.mean);
perfMean = mean(Qperf);
perfSEM = std(Qperf)/sqrt(size(Qperf,1));

Qrt = cat(1,RT.mean);
rtMean = mean(Qrt);
rtSEM = std(Qrt)/sqrt(size(Qrt,1));


ntrial_bin = edges(end)/nbins(1); % number of trials per bin


%% Plot performance
width = 16;
height = 12;
LW = 1.5;
FS = 16;
figure(1)
clf;

xBnd = [0.5 sum(nbins)+0.5];
yBnd = [0.08 .15];

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

hold on
%fill([xBnd(1) nbins(1)+0.5 nbins(1)+0.5 xBnd(1)], [yBnd(1) yBnd(1) yBnd(2) yBnd(2)],'k','FaceColor',[1 1 1]*1);
fill(ntrial_bin*[0.5 nbins(1)+0.5 nbins(1)+0.5 0.5], [yBnd(1) yBnd(1) yBnd(2) yBnd(2)],...
    'k','FaceColor',[1 1 1]*0.9,'EdgeColor','none');
fill(ntrial_bin*[nbins(1)+0.5 sum(nbins)+0.5 sum(nbins)+0.5 nbins(1)+0.5], [yBnd(1) yBnd(1) yBnd(2) yBnd(2)],...
    'k','FaceColor',[1 1 1]*0.8,'EdgeColor','none');
%plot(binCenter,Qdev','Color',[1 1 1]*0.5);
errorbar(ntrial_bin*(1:length(perfMean)),perfMean,perfSEM,'Color','k','LineWidth',LW,'LineStyle','-');

xlim(xBnd*ntrial_bin);
ylim(yBnd);
xlabel('trial', 'FontSize', FS, 'FontName', 'Times');
ylabel('absolute deviation from optimal', 'FontSize', FS, 'FontName', 'Times');

% lh = legend(h, legLabel);
% set(lh, 'Location', 'northwest', 'FontSize', 12, 'FontName', 'Times');

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');

% Inset
inset = axes('Position',[0.52 0.27 0.35 0.2]);
hold on
SEM = std(perfFirst,'omitnan')/sqrt(numel(subInd));
%errorbar(1:nfirst,mean(devFirst,'omitnan'),std(devFirst,'omitnan')/sqrt(numel(subInd)),'Color','k','LineWidth',LW,'LineStyle','-');
plot(1:nfirst,mean(perfFirst,'omitnan'),'Color','k','LineWidth',LW,'LineStyle','-');
plot(1:nfirst,mean(perfFirst,'omitnan')+SEM,'Color','r','LineWidth',1,'LineStyle','-.');
plot(1:nfirst,mean(perfFirst,'omitnan')-SEM,'Color','r','LineWidth',1,'LineStyle','-.');

xlim([1 nfirst]);
ylim([0.05 .22]);

set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times','XTick',[1 25 50]);

% Print
print(gcf, '-dpng', '-r400', '.\..\..\..\plots\exp2\supp_info\performance_trials_prior.png');
