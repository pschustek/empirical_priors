% Fig. S3
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = [1:20 22:24];

% Load data
load('.\..\..\..\models\exp2\model_evidence')

M = M(subInd,:);
rowNames = {rowNames{subInd}};

%tM = array2table(M,'RowNames',rowNames,'VariableNames',fnames)

%% Select (order of) models 
list = {'opt_zmap','tly_zmap','avg_zmap','diff_zmap'};
modelLabel = {'tly','avg','diff'};

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
xrng = [.3 3.7];

% Plots
set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

hold on
bar(expected,'FaceColor',hsv2rgb([0 .5 0.9]),'LineWidth',1.2,'BarWidth',0.5);
for j=1:numModels-1
    line(j*[1 1],expectedCI(:,j),'LineWidth',2,'Color','k');
    if round(exceed,3) < 1
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
print(gcf, '-dpng', '-r400', '.\..\..\..\plots\exp2\si_prior_bms_paired_ref_model_comparison.png');