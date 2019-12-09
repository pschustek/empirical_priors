% Supp Fig. 1b
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = 1:24;

% Load data
load('.\..\..\..\data\exp2_data.mat');

% load fitted models
load('.\..\..\..\models\exp2\prior_cvll.mat');

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

%%
for s=subInd
    clear trials
    trials = trialData{s};
    
    C = trials.confDec;
    Y = trials.confHeads;
    mEv = trials.meanEvidence;
    N = trials.sampleSize; nH = trials.meanEvidence.*N;
    
    % Objective probability of being correct for decision [1,0]
    objDecCorr = [trials.optConfHeads 1-trials.optConfHeads];
    
    % Probability correct for participant's response
    probDecCorr = sum(objDecCorr.*[trials.decision ~trials.decision],2);
    
    optC = trials.optConfHeads; % response of the optimal model
    optCA = abs(optC-0.5) + 0.5; % aligned with larger proportion
    Y(optC<0.5) = 1-Y(optC<0.5); % aligned participant response with larger proportion
    optCorr = trials.optDecHeads == trials.coinBias;
    numTrials = size(trials.num,1);
    %  fitted optimal model with distortion
    distC = cvll_prior(s).opt_opt_sig.Y;
    distC(optC<0.5) = 1-distC(optC<0.5);
    distCorr = (cvll_prior(s).opt_opt_sig.Y>.5) == trials.coinBias;
    
    
    % Optimal model
    binIdx = discretize(optCA,edges); % aligned confidence
    mask = binIdx==[1:nbins];
    for j=1:size(mask,2)
        opt(s).resp_vs_opt(j) = mean(optCA(mask(:,j)),'omitnan'); % actual response vs optimal response
        sub(s).resp_vs_opt(j)  = mean(Y(mask(:,j)),'omitnan'); % actual response vs optimal response
        dist(s).resp_vs_opt(j)  = mean(distC(mask(:,j)),'omitnan'); % distorted response vs optimal response
    end
    
    
end

%% Aggregate across participants


Map_opt = cat(1,opt.resp_vs_opt);
Map = cat(1,sub.resp_vs_opt);
Map_dist = cat(1,dist.resp_vs_opt);


Mapmean = mean(Map,'omitnan');
MapSEM = std(Map,'omitnan')/sqrt(size(Map,1));
Map_distmean = mean(Map_dist,'omitnan');


% Tests of deviation from unity line
for j=1:size(Map,2)
    tmp = Map(:,j);
    pval.above(j) = signrank(tmp(~isnan(tmp)),binCenter(j),'tail','right');
    pval.below(j) = signrank(tmp(~isnan(tmp)),binCenter(j),'tail','left');
end

% Group level calibratedness: Correlation with model
mask = ~isnan(Map);
[rho,pcal] = corr(Map(mask),Map_opt(mask));
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

q = 0.5:0.01:1;
plot(q,q,'Color',[1 1 1]*0.6,'LineWidth',0.8,'LineStyle','--');
h(2) = errorbar(binCenter, Mapmean',MapSEM','Color','k','LineWidth',LW,'LineStyle','-','CapSize',3);
h(3) = plot(binCenter,Map_distmean,'Color','g');

for j=1:4
    if pval.above(j)<=0.05 && pval.above(j)>0.01
        text(binCenter(j),Mapmean(j)+0.07,'$\ast$','Interpreter','latex','FontSize', FS-1,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    if pval.above(j)<=0.01
        text(binCenter(j),Mapmean(j)+0.07,'$\ast\ast$','Interpreter','latex','FontSize', FS-1,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

for j=5:length(binCenter)
    if pval.below(j)<=0.05 && pval.below(j)>0.01
        text(binCenter(j),Mapmean(j)-0.07,'$\ast$','Interpreter','latex','FontSize', FS-1,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    if pval.below(j)<=0.01
        text(binCenter(j),Mapmean(j)-0.07,'$\ast\ast$','Interpreter','latex','FontSize', FS-1,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

xlim([0.5 1]);
ylim([0.5 1]);

xlabel('optimal response', 'FontSize', FS, 'FontName', 'Times');
ylabel('average subject response', 'FontSize', FS, 'FontName', 'Times');

legend(h(2:3), {'subject', 'optimal with distortion'}, 'location', 'northwest');

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','YTick',0.5:0.1:1,'XTick',0.5:0.1:1);


%% Print
print(gcf, '-dpng', '-r400', '.\..\..\..\plots\exp2\supp_info\si_prior_confidence_calibration2.png');



