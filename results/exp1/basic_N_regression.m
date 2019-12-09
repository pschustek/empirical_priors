% Fig. 2b, Supp Fig 2 and Supp Fig 3
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\src\');

% Subselect subjects
subInd = 1:24;

% Load data
load('.\..\..\data\exp1_data.mat');

% Load model fits
load('.\..\..\models\exp1\basic_cvll.mat');

modelnames = {'mEv','diff'}; % models: ratio and difference
n_models = length(modelnames);

nPerm = 10000;

%% Compute pattern
% Takes ~ 10 min
tic;
for s=subInd % for each subject
    
    clear trials
    trials = trialData{s};
    fprintf('- computing pattern for participant %d\n', s);
    
    mEv = trials.meanEvidence;
    N = trials.sampleSize;
    numTrials = size(trials.num,1);
    
    response = trials.confHeads; % subject response
    optRes = trials.optConfHeads; % optimal model response
    
    %% Fit N
    setN = unique(N); % sample sizes used
    x0 = 0; 5;
    lb = -1E2;
    ub = 1E2;
    opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    
    % Normalize mEv
    mEv = 2*(mEv-0.5);
    
    slope = nan(length(setN),1);
    R2exp = nan(length(setN),1);
    slope_opt = nan(length(setN),1);
    R2opt = nan(length(setN),1);
    R2mod = nan(length(setN),n_models);
    
    for j=1:length(setN) % for each sample size
        mask = setN(j)==N; % select trials with given sample size
        ff = @(w,X) genlogf(X,w); % logistic function (only parameter is slope)
        
        % Participant
        fObj = @(w) -f_obj( ff(w,mEv(mask)), response(mask) );
        [slope(j),~,exitflag] = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
        R2exp(j) = rsquared( ff(slope(j), mEv(mask)), response(mask) );
        
        % Model
        fObj = @(w) -f_obj( ff(w,mEv(mask)), optRes(mask) );
        [slope_opt(j),~,exitflag] = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
        R2opt(j) = rsquared( ff(slope_opt(j), mEv(mask)), optRes(mask) );
                   
        % slope for fitted models
        for m=setdiff(1:n_models,5)
            strk =  cvll_basic(s).(modelnames{m}); % model
            fObj = @(w) -f_obj( ff(w,mEv(mask)), strk.Y(mask) );
            [slope_mod(j,m),~,exitflag] = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
            R2mod(j,m) = rsquared( ff(slope(j), mEv(mask)), response(mask) );
        end
        
    end
    
    sub(s).slope = slope;
    sub(s).slope_opt = slope_opt;
    sub(s).R2exp = R2exp;
    sub(s).R2opt = R2opt;
    sub(s).slope_mod = slope_mod;
    sub(s).R2mod = R2mod;
    
    % Individually increasing slope with N
    linMod = fitlm(setN,slope,'linear');
    sub(s).rho_lin_N = linMod.Coefficients.Estimate(2);
    sub(s).p_lin_N = linMod.coefTest;
    
    % Permutation test; Shuffle N
    % Null-hypothesis: N does not have an effect
    res = nan(nPerm,1);
    parfor j=1:nPerm
        idx = randperm(size(slope,1));
        
        mdl = fitlm(setN(idx), slope);
        
        res(j) = mdl.Coefficients.Estimate(2);
    end
    sub(s).p_lin_N_perm = mean(res>=sub(s).rho_lin_N);
    
end
toc;
%% Aggregate across participants
allw = cat(2,sub.slope);
allOptw = cat(2,sub.slope_opt);

wMean = mean(allw,2);
wSEM = std(allw,[],2)/sqrt(numel(subInd));

sMeanOpt = mean(cat(2,sub.slope_opt),2);
sMeanMod = mean(cat(3,sub.slope_mod),3);


%% Group level (regression) to evidence positive slope
tmp = repmat(setN,1,numel(subInd));
[rho, p] = corr(tmp(:),allw(:));
fprintf('- [result] correlation N with slope, r = %.3f, p = %.3e\n', rho, p);

%% Number of participants for which relationship is individually significant
numSignificant = sum(cat(1,sub.p_lin_N_perm)<0.05);
fprintf('- [result] individually significant N-dependence = %d out of %d\n', numSignificant, numel(subInd));

%% test deviation from optimal model at N=13
[~, pval_N13] = ttest(allw(end,:), allOptw(end,:));
fprintf('- [result] T-test of participants vs optimal model at sample size=13: p = %f\n', pval_N13);


%% Plot
for i=1:n_models+1
    figure(i); 
    width = 8;
    height = 6;
    LW = 1.2;
    FS = 11;
    clf;
    
    hold on
    
    colGrad = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*.5, ones(numel(q),1)*0.85]);
    cmap = colGrad(0:0.02:1);
    colormap(cmap);
    if i==1
        patch([setN' nan],[sMeanOpt' nan],[linspace(0,1,numel(setN)) nan],'EdgeColor','interp','Marker','none','LineWidth', LW+2, 'linestyle','-');
     %   patch([setN' nan],[sMeanHeur' nan],[linspace(0,1,numel(setN)) nan],'EdgeColor','interp','Marker','none','LineWidth', LW+2);
        set(gcf, 'name','regression exp1');
    else
        patch([setN' nan],[sMeanMod(:,i-1)' nan],[linspace(0,1,numel(setN)) nan],'EdgeColor','interp','Marker','none','LineWidth', LW+2);
        set(gcf, 'name',sprintf('regression exp1 %s',modelnames{i-1}));
    end
    % h(1) = plot(setN,sMeanOpt,'Color',[1 0.5 0.5], 'LineWidth', LW+2);
    h(2) = errorbar(setN,wMean,wSEM,'Color','k','LineWidth',LW,'LineStyle','none');
    
    % Try
    colGrad = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*1, ones(numel(q),1)*0.85]);
    cmap = colGrad(0:0.02:1);
    colormap(cmap);
    patch([setN' nan],[wMean' nan],[linspace(0,1,numel(setN)) nan],'EdgeColor','interp','Marker','o',...
        'MarkerEdgeColor', 'k', 'LineWidth', .9,'MarkerFaceColor','flat','LineStyle', 'none');
    
    xlim([2 14]);
    ylim([0 10]);
    xlabel('sample size','Interpreter','latex');
    ylabel('slope confidence curve');
    
    % lh = legend(h,'optimal','experiment','Location','best');
    % set(lh,'Interpreter','latex','Box','off')
    
    set(gcf,'Color',[1,1,1]);
    
    % Position plot on the screen for drawing
    set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);
    
    % Position plot on the paper for printing
    set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
        'PaperSize', [width height], 'PaperPosition', [0 0 width height]);
    
    % Axes
    set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
        'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',3:2:13);
    
   % figure name
    figname = 'basic_N_regression';
    if i>1
        figname = [figname '_' modelnames{i-1}];
    end
    set(gcf,'name',figname);
    
    
    %% Print
    filename =  fullfile('.\..\..\plots\exp1\',[figname '.png']);
    print(gcf, '-dpng', '-r400', filename);
end  
    

%% plot individual response only
    figure(2*n_models+2); 
    clf;
    LW = .8;
    
    for s=subInd
        subplot(4,6,s); title(sprintf('S%d',s)); hold on;
        set(gca, 'fontsize',5);
        
        % Try
        colGrad = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*1, ones(numel(q),1)*0.85]);
        cmap = colGrad(0:0.02:1);
        colormap(cmap);
        patch([setN' nan],[allw(:,s)' nan],[linspace(0,1,numel(setN)) nan],'EdgeColor','interp','Marker','o',...
            'MarkerEdgeColor', 'k', 'LineWidth', .9,'MarkerFaceColor','flat','LineStyle', 'none', 'markersize',3);
        
        if s==22
            xlabel('sample size','Interpreter','latex', 'fontsize',9);
        elseif s<19
            set(gca, 'xticklabel','');
        end
        if s==13
            ylabel('slope confidence curve', 'fontsize',9);
        end
        axis tight;
    end    
    
    set(gcf,'Color',[1,1,1]);
    
    % Position plot on the paper for printing
    set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
        'PaperSize', [width height], 'PaperPosition', [0 0 width height]);
    
    % save figure
    figname = 'basic_N_regression_individual_subj';
    set(gcf,'name',figname);
    
    % Print
    filename =  fullfile('.\..\..\plots\exp1\',[figname '.png']);
    print(gcf, '-dpng', '-r400', filename);


       filename =  fullfile('.\..\..\plots\exp1\',[figname '.pdf']);
    print(gcf, '-dpdf', '-r400', filename);