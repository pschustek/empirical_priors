% Supp Fig. 4
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = 1:24;


% Load data
load('.\..\..\..\models\exp1\model_evidence');

% CVLL for models and individual participants
M = M(subInd,:);
rowNames = {rowNames{subInd}};
tM = array2table(M,'RowNames',rowNames,'VariableNames',fnames)

%% Select (order of) models
list = {'opt','mEv','diff'};

modelLabel = {'ratio','difference'};

IDX = nan(length(list),1);
for m=1:length(list)
    IDX(m) = find(strcmp(fnames, list{m}));
end

%% Pairwise BMS comparison
numIter = 50;
numModels = length(IDX);

medDiff = nan(1,numModels-1);
exceed = nan(1,numModels-1);
expected = nan(1,numModels-1);
expectedCI = nan(2,numModels-1);

for j=2:numModels
    comp = [IDX(1) IDX(j)];         % compare only pairs
    alpha = bms(numIter,M(:,comp));
    
    fh = @(q) bms(numIter,q);
    
    % Bootstrapped BMS
    bs = bootstrp(1000,fh,M(:,comp));
    R = bs./repmat(sum(bs,2),1,2);
    
    % Expected likelihood of model k
    r = alpha./sum(alpha);
    expected(1,j-1) = r(1);
    expectedCI(:,j-1) = quantile(R(:,1),[0.05 0.95])';
    
    % Standard difference
    medDiff(1,j-1) = -median(dHart(diff(M(:,comp),1,2)));
    
    % Beta-marginal
    exceed(1,j-1) = betacdf(0.5,alpha(1),sum(alpha)-alpha(1),'upper');       % actual exceedance probability: 1 vs. the other
end

%% Plot: Expectation Values
figure(1);
clf;

width = 7.5;
height = width/1.25;
FS = 9;
LW = 1.2;
xrng = [.3 2.7];

% Plots
set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

hold on
bar(expected,'FaceColor',hsv2rgb([0 .5 0.9]),'LineWidth',1.2,'BarWidth',0.5);
% errorbar(1:numModels-1,expected,expected-expectedCI(1,:),expected-expectedCI(2,:),'LineStyle','none','CapSize',0,'Color','k','LineWidth',2);
for j=1:numModels-1
    line(j*[1 1],expectedCI(:,j),'LineWidth',2,'Color','k');
    if round(exceed(j),3) < 1
        text(j,max(expectedCI(2,:))+0.06,['$p_{e}>' num2str(round(exceed(j),3)) '$'],'Interpreter','latex','FontSize', FS,...
            'FontName','Times','HorizontalAlignment','center','HorizontalAlignment','center');
    else
        text(j,max(expectedCI(2,:))+0.06,['$p_{e}\approx 1$'],'Interpreter','latex','FontSize', FS,...
            'FontName','Times','HorizontalAlignment','center','HorizontalAlignment','center');
    end
end

xlim(xrng);
ylim([.5 1])
ylabel('$p(\mathrm{opt})$','FontSize', FS,'FontName','Times','Interpreter','latex');

set(gca,'XTick',1:numModels-1,'YTick',0.5:0.1:1,...
    'FontSize', FS, 'FontName', 'Times','Color','none','TickLabelInterpreter','latex',...
    'Position',[0.2 0.12 0.78 0.75],'XTickLabel',modelLabel);


%% Print
print(gcf, '-dpng', '-r400', '.\..\..\..\plots\exp1\supp_info\si_basic_bms_paired_ref_model_comparison.png');


%% plot histogram of difference
Models_incl = 1:2;
figure(2);
width = 2*8;
height = 6.5;
LW = 1.2;
FS = 11;
clf;

M_diff = M(:,IDX(1)) - M(:,IDX(2:end));
binsize = 5;
bins = binsize*(floor(min(M_diff(:))/binsize):ceil(max(M_diff(:))/binsize));
for j=1:length(Models_incl)
    jj = Models_incl(j);
    subplot(1,length(Models_incl),j); hold on; title(strrep(modelLabel{jj},'_',''));
    xx = M_diff(:,jj);
    histogram(xx(xx<=-5),bins, 'facecolor','r');
    histogram(xx(xx>-5 & xx<5),bins, 'facecolor',.7*[1 1 1]);
    histogram(xx(xx>=5),bins, 'facecolor','b');
    
    ylim([0 14.2]);
    plot([0 0], ylim,'r','linewidth',2);
    if j==1
        ylabel('number of participants');
        %  xlabel('Difference in model evidence \Delta');
        text(20, -3, 'Difference in model evidence \Delta (dHart)'); % custom made xlabel
    end
    set(gca, 'Position', get(gca, 'Position')+[0 .1 0 -.1]);
    
    % Axes
    set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
        'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');
end

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

%% Print
print(gcf, '-dpng', '-r400', '.\..\..\..\plots\exp1\supp_info\si_basic_bms_paired_ref_model_comparison_individual.png');