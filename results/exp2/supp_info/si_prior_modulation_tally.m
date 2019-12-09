% Supp Fig 10 a.e

clear variables
close all

nu1 = 12;
nu2 = 9; % for heuristics model

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = [1:24]; 

% Load data
load('.\..\..\..\data\exp2_data.mat');

% Load model fits
load('.\..\..\..\models\exp2\prior_cvll.mat');
modelnames = {'avg_zmap','diff_zmap','tly_zmap'};
n_models = length(modelnames);


% Make fixed mEv bins
nbins = 4;
edges = [0 .45 .55 .65 .72 1];
binCenter = [diff(edges)/2 + edges(1:end-1)]';

%% Compute pattern
for s=subInd
    
    clear trials
    trials = trialData{s};
    
    % compute tally
    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    numTrials = length(trials.num);
    block = ceil([1:numTrials]'/unique(trials.blockLength));
    tally =  estimate_M.tally( [mEv, N, block]);
    
    
    response = trials.confHeads;
    optRes = trials.optConfHeads;
    
    inBlock = trials.trialInBlock;
    
  
    % Align with real block bias
    blockbias = 2*(trials.blockBias-0.5);
    response = blockbias.*(response-0.5) + 0.5;
    optRes = blockbias.*(optRes-0.5) + 0.5;
    mEv = blockbias.*(mEv-0.5) + 0.5;
    tally = blockbias.*(tally-0.5) + 0.5;
    
    % Remove first in-block trial
    idx = (inBlock==1); %1:trials.blockLength:numel(mEv);
    response(idx) = [];
    optRes(idx) = [];
    mEv(idx) = [];
    tally(idx) = [];
    inBlock(idx) = [];
    
    %% Make groupings according to tally
    % Mask for mEv
    binIdx = discretize(tally,edges);
    
      
    [G,bb, bc] = findgroups(binIdx, inBlock);
    res = splitapply(@mean,response,G);
    sub(s).confBin = reshape(res, 4, length(edges)-1);
    
    res = splitapply(@mean,optRes,G);
    opt(s).confBin = reshape(res, 4, length(edges)-1);
    
    % fitted models
    for m=1:n_models
        strk =  cvll_prior(s).(modelnames{m}); % model
        modRes = blockbias.*(strk.Y-0.5) + 0.5;
        modRes(idx) = [];
        res = splitapply(@mean,modRes,G);
        modl(s,m).confBin = reshape(res, 4, length(edges)-1);
    end 
    
    
end

%% Aggregate across participants
confMean = mean(cat(3,sub.confBin),3);
confSEM = std(cat(3,sub.confBin),[],3)/sqrt(numel(subInd));

confMeanOpt = mean(cat(3,opt.confBin),3);

for m=1:n_models
    confMeanMod{m} = mean(cat(3,modl(:,m).confBin),3);
end



%% Plot
for i=1:n_models+2
    figure(i);
    m= i-1;
    figname = 'prior_modulation_tally';
    if i==n_models+2
        figname = [figname '_subject'];
    elseif i>1
        figname = [figname '_' modelnames{m}];
    end
    set(gcf, 'name',figname);
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
    
    colormap = cool;
    colors = colormap([4 24 44 64],:);
    
    % Model
    for j=1:4
        blockLineStyle = {'-','-','-'};
        colorMap = [hsv2rgb([.08 .5 0.85]); hsv2rgb([.33 .5 .85]) ; 0.59 0.82 0.59];
        
        if i==1
            h(j,1) = plot(binCenter*100,confMeanOpt(j,:),'Color',colors(j,:),'LineStyle',blockLineStyle{1}, 'LineWidth', LW+2);
        elseif i<n_models+2
            h(j,1) = plot(binCenter*100,confMeanMod{m}(j,:),'Color',colors(j,:),'LineStyle',blockLineStyle{1}, 'LineWidth', LW+2);
        else
            colorMap = [hsv2rgb([.08 1 0.85]); hsv2rgb([.33 1 .85]); nan(1,3)];
            h(j,1) =    errorbar(binCenter*100,confMean(j,:),confSEM(j,:),'Color',colors(j,:),'LineStyle','-', 'LineWidth', LW+2);
        end
        
        % Participants
          leg{j} = ['trial position ' num2str(j+1)];
    end
    
   
    xlim([0.2 .9]*100);
    ylim([0.4 0.82]);
    
    clear xL
    xL = xlabel('aligned cumulative proportion','interpreter','latex','FontSize',11);
    % set(xL,'Position',[53 xL.Position(2:end)]);
    ylabel('aligned confidence');
    
    if i==1
        legh = legend(h, leg,'location','best');
        set(legh, 'location','none','position',[.55 .32 .38 .16]);
    end
      
    set(gcf,'Color',[1,1,1]);
    
    % Position plot on the screen for drawing
    set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);
    
    % Position plot on the paper for printing
    set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
        'PaperSize', [width height], 'PaperPosition', [0 0 width height]);
    
    % Axes
    set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...
        'XMinorTick', 'off', 'YMinorTick', 'off', 'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');
    
    %% save figure
    filename = fullfile('.\..\..\..\plots\exp2\supp_info\', [figname '.png']);
    print(gcf, '-dpng', '-r400', filename);
    
    
end


