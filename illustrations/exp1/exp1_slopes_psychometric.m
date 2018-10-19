% Fig. 1d
clear variables
close all

%% Prepare
addpath('.\..\..\src\');

numTrials = 50000;
numDist.support = 3:13;
numMap = @(n) 0*n + 1;
numDist.pmf = numMap(numDist.support)/sum(numMap(numDist.support));
tmp = repmat(numDist.support,numTrials,1);
N = tmp(logical(mnrnd(1,numDist.pmf,numTrials)));
mu = betarnd(4,4,numTrials,1);
nH = binornd(N,mu);
mEv = nH./N;

optRes = opt_inf.basic_confH(nH,N,4,4);


%% Fit N
setN = unique(N);
x0 = [5];
lb = -1E2*[1];
ub = 1E2*[1];
opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');

% Normalize mEv
mEv = 2*(mEv-0.5);

slope_opt = nan(length(setN),1);
R2opt = nan(length(setN),1);
for j=1:length(setN)
    mask = setN(j)==N;
    ff = @(w,X) genlogf(X,w);
    
    % Model
    fObj = @(w) -f_obj( ff(w,mEv(mask)), optRes(mask) );
    [slope_opt(j),~,exitflag] = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
    R2opt(j) = rsquared( ff(slope_opt(j), mEv(mask)), optRes(mask) );
    
end

%% Plot
figure(1);
width = 8;
height = 5.5;
LW = 1.5;
FS = 11;
clf;

colGrad = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*1, ones(numel(q),1)*0.85]);
cmap = colGrad(0:0.02:1);
colormap(cmap);

hold on

% plot(setN,slope_opt,'Color','k', 'LineWidth', LW);
patch([setN' nan],[slope_opt' nan],[linspace(0,1,numel(setN)) nan],'EdgeColor','interp','Marker','none','MarkerFaceColor','flat','LineWidth', LW+2);

xlim([2 14]);
ylim([0 5.5]);
xlabel('sample size','Interpreter','latex');
ylabel('slope confidence curve','Interpreter','latex');

% lh = legend(h,'optimal','experiment','prior mismatch','Location','best');
% set(lh,'Interpreter','latex','Box','off')

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',3:2:13, 'YTick', 0:2.5:5);

%% Print
print(gcf, '-dpng', '-r400', 'N_slope_psychometric.png');