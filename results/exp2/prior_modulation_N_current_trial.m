% Fig. 6b
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\src\');

% Subselect subjects
subInd = [1:20 22:24];

% Load data
load('.\..\..\data\exp2_data.mat')

% Make fixed mEv bins
nbins = 5;
% Pool all responses together
accum = [];
for s=subInd
    mEv = trialData{s}.meanEvidence;
    bblock = 2*(trialData{s}.coinBias-0.5);
    accum = [accum; bblock.*(mEv-0.5) + 0.5];
end
% With respect to equally filled bins
qtls = [0:nbins]/nbins;
edges = quantile(accum,qtls(2:end-1));
edges = [0 edges 1];
binCenter = [diff(edges)/2 + edges(1:end-1)]';

%% Prepare data
for s=subInd

    clear trials    
    trials = trialData{s};

    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    
    response = trials.confHeads;
    optRes = trials.optConfHeads;

    % Align with real block bias
    block = 2*(trials.blockBias-0.5);
    response = block.*(response-0.5) + 0.5;
    optRes = block.*(optRes-0.5) + 0.5;
    mEv = block.*(mEv-0.5) + 0.5;

    %% Make groupings according to mEv and previous block evidence
    % Mask for mEv 
    binIdx = discretize(mEv,edges);     

    few = N < quantile(N,0.4);              
    many = N >= quantile(N,0.6);
    groupMask = logical([few many]);    
    allConditions = any(groupMask,2);   
    
    [G,~,bb] = findgroups(binIdx(allConditions), groupMask(allConditions,1));
    res = splitapply(@mean,response(allConditions),G);
    sub(s).confBin = [res(bb==1) res(bb==0)];
    
    res = splitapply(@mean,optRes(allConditions),G);
    opt(s).confBin = [res(bb==1) res(bb==0)];
    
    % Fit for small N
    mdl = fitlm(mEv(few), response(few));
    sub(s).few_anova_pval = mdl.anova.pValue(1);
    sub(s).few_slope = mdl.Coefficients.Estimate(2);
    
    % Fit for large N
    mdl = fitlm(mEv(many), response(many));
    sub(s).many_anova_pval = mdl.anova.pValue(1);
    sub(s).many_slope = mdl.Coefficients.Estimate(2);
    
    sub(s).diff_slope = sub(s).many_slope - sub(s).few_slope;

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

% Crossover pattern (slope difference > 0)
[p,~,stats] = signrank(vertcat(sub.diff_slope),0,'tail','right');
fprintf('- [result] p-value = %.3e\n', p);

%% Plot
figure(1);
clf;
LW = 1.2;
FS = 7;
width = 5.2;
height = width/1.333;
clf;
hold on

clear h hl

line([0 100],0.5*[1 1],'LineStyle','--','LineWidth',0.8,'Color',[1 1 1]*0.8);
line(50*[1 1],[0 1],'LineStyle','--','LineWidth',0.8,'Color',[1 1 1]*0.8);

% Optimal model
blockLineStyle = {'-','-','-'};
colorMap = [hsv2rgb([.08 .5 0.85]); hsv2rgb([.33 .5 .85]) ; 0.59 0.82 0.59];
legIdx = 1:2;
for b=1:size(groupMask,2)
    h(legIdx(b)) = plot(binCenter*100,confMeanOpt(:,b),'Color',colorMap(b,:),'LineStyle',blockLineStyle{b},'LineWidth',LW+3);
end

% Participants
blockLineStyle = {'-','-','-'};
colorMap = [hsv2rgb([.08 1 0.85]); hsv2rgb([.33 1 .85]); nan(1,3)];
legIdx = 1:2;
for b=1:size(groupMask,2)
    errorbar(binCenter*100,confMean(:,b),confSEM(:,b),'Color',colorMap(b,:),'LineStyle',blockLineStyle{b},'LineWidth',LW);
end

for j=1:1
    if all_pval_few(j)<=0.05 & all_pval_few(j)>0.01
        text(binCenter(j)*100,confMean(j)+0.05,'$\ast$','Interpreter','latex','FontSize', FS,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    if all_pval_few(j)<=0.01 
        text(binCenter(j)*100,confMean(j)+0.05,'$\ast\ast$','Interpreter','latex','FontSize', FS,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

for j=3:length(binCenter)
    if all_pval_many(j)<=0.05 & all_pval_many(j)>0.01
        text(binCenter(j)*100,confMean(j)-0.05,'$\ast$','Interpreter','latex','FontSize', FS,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    if all_pval_many(j)<=0.01 
        text(binCenter(j)*100,confMean(j)-0.05,'$\ast\ast$','Interpreter','latex','FontSize', FS,'FontName','Times',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

xlim([15 90]);
ylim([0.25 1]);

xlabel('aligned sample proportion (\%)','interpreter','latex');
ylabel('aligned confidence');

% Custom legend
xl = xlim;
yl = ylim;
xpos = 20;
ypos = 0.8;
rectangle('Position',[xpos ypos diff(xl)/10 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.33 1 .85]));
ypos = 0.7;
rectangle('Position',[xpos ypos diff(xl)/10 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.08 1 0.85]));

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XMinorTick', 'off', 'YMinorTick', 'off', 'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');

%%
print(gcf, '-dpng', '-r400', '.\..\..\plots\exp2\prior_modulation_N_current_trial.png');