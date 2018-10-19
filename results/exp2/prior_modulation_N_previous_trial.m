% Fig. 6c
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\src\');

% Subselect subjects
subInd = [1:20 22:24];

% Load data
load('.\..\..\data\exp2_data.mat');

% Make fixed mEv bins
nbins = 4;
edges = [0 0.4545 0.6250 0.7500 1.0000]; 
binCenter = [diff(edges)/2 + edges(1:end-1)]';

%% Compute pattern
for s=subInd

    clear trials    
    trials = trialData{s};

    mEv = trials.meanEvidence;

    % Block evidence prior to present trial
    pmEv = [nan; mEv(1:end-1)];
    pmEv(1:trials.blockLength:numel(mEv)) = nan;
    pN = [nan; trials.sampleSize(1:end-1)];
    pN(1:trials.blockLength:numel(mEv)) = nan;

    response = trials.confHeads;
    optRes = trials.optConfHeads;

    % Align with real block bias
    block = 2*(trials.blockBias-0.5);
    response = block.*(response-0.5) + 0.5;
    optRes = block.*(optRes-0.5) + 0.5;
    mEv = block.*(mEv-0.5) + 0.5;
    pmEv = block.*(pmEv-0.5) + 0.5;
    
    % Remove first in-block trial
    idx = 1:trials.blockLength:numel(mEv);
    response(idx) = [];
    optRes(idx) = [];
    mEv(idx) = [];
    pmEv(idx) = [];
    pN(idx) = [];

    %% Make groupings according to previous sample size
    % Mask for mEv 
    binIdx = discretize(pmEv,edges);
    
    few = pN <= quantile(pN,0.5);
    many = pN > quantile(pN,0.5);
    groupMask = logical([few many]);    
    allConditions = any(groupMask,2);   
    
    [G,~,bb] = findgroups(binIdx(allConditions), groupMask(allConditions,1));
    res = splitapply(@mean,response(allConditions),G);
    sub(s).confBin = [res(bb==1) res(bb==0)];
    
    res = splitapply(@mean,optRes(allConditions),G);
    opt(s).confBin = [res(bb==1) res(bb==0)];
    
    % Higher modulation for higher sample size
    X = pmEv;
    X_F = X(few);
    X_M = X(many);
    Y_F = response(few);
    Y_M = response(many);
    mdl_F = fitlm(X_F, Y_F);
    mdl_M = fitlm(X_M, Y_M);
    sub(s).slope_few = mdl_F.Coefficients.Estimate(2);
    sub(s).slope_many = mdl_M.Coefficients.Estimate(2);

end

%% Aggregate across participants
confMean = mean(cat(3,sub.confBin),3);
confSEM = std(cat(3,sub.confBin),[],3)/sqrt(numel(subInd));

confMeanOpt = mean(cat(3,opt.confBin),3);

% Test on group level
tmp = cat(1,sub.confBin);
tmpF = reshape(tmp(:,1),nbins,[]);
tmpM = reshape(tmp(:,2),nbins,[]);
for n=1:nbins
    F = tmpF(n,:);
    M = tmpM(n,:);
    % many > few
    all_pval_many(n) = signrank(M,F,'tail','right');
    % few > many
    all_pval_few(n) = signrank(M,F,'tail','left');
end

% Higher modulation for higher sample size (group level)
[p,~,stats] = signrank(vertcat(sub.slope_many),vertcat(sub.slope_few),'tail','right');
fprintf('- [result] p-value = %.4f\n', p);

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

% Model
blockLineStyle = {'-','-','-'};
colorMap = [hsv2rgb([.08 .5 0.85]); hsv2rgb([.33 .5 .85]) ; 0.59 0.82 0.59];
legIdx = 1:2;
for b=1:size(groupMask,2)
    h(legIdx(b)) = plot(binCenter*100,confMeanOpt(:,b),'Color',colorMap(b,:),'LineStyle',blockLineStyle{b}, 'LineWidth', LW+2);
end

% Participants
blockLineStyle = {'-','-','-'};
colorMap = [hsv2rgb([.08 1 0.85]); hsv2rgb([.33 1 .85]); nan(1,3)];
legIdx = 1:2;
for b=1:size(groupMask,2)
    errorbar(binCenter*100,confMean(:,b),confSEM(:,b),'Color',colorMap(b,:),'LineStyle',blockLineStyle{b}, 'LineWidth', LW);
end

for j=1:2
    if all_pval_few(j)<=0.05 & all_pval_few(j)>0.01
        text(binCenter(j)*100,confMean(j)+0.03,'$\ast$','Interpreter','latex','FontSize', FS,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    if all_pval_few(j)<=0.01 
        text(binCenter(j)*100,confMean(j)+0.03,'$\ast\ast$','Interpreter','latex','FontSize', FS,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

for j=4:length(binCenter)
    if all_pval_many(j)<=0.05 & all_pval_many(j)>0.01
        text(binCenter(j)*100,confMean(j)-0.03,'$\ast$','Interpreter','latex','FontSize', FS,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    if all_pval_many(j)<=0.01 
        text(binCenter(j)*100,confMean(j)-0.03,'$\ast\ast$','Interpreter','latex','FontSize', FS,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

xlim([0.2 .9]*100);
ylim([0.55 0.82]);

clear xL
xL = xlabel('aligned proportion previous sample (\%)  ','interpreter','latex','FontSize',11);
% set(xL,'Position',[53 xL.Position(2:end)]);
ylabel('aligned confidence');

% Custom legend
xl = xlim;
yl = ylim;
xpos = 24;
ypos = 0.76;
rectangle('Position',[xpos ypos diff(xl)/10 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.33 1 .85]));
ypos = 0.72;
rectangle('Position',[xpos ypos diff(xl)/10 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.08 1 0.85]));

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
print(gcf, '-dpng', '-r400', '.\..\..\plots\exp2\prior_modulation_N_previous_trial.png');