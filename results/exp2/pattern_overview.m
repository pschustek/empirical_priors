%% Illustrate all patterns for optimal inference 
% Fig. 4

clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\src\');

rng(87,'twister');

%% Simulate behavior
blockLength = 5;
numSim = 1000000;
[ confHeads, mEv, N, blockBias, pbEv, bEv, coinBias ] = sim_opt_inf(numSim);

% Important quantities
inBlock = mod([1:length(confHeads)]'-1,5) + 1;
blockIdx = ceil([1:numSim]/blockLength)';
confDec = abs(confHeads-0.5) + 0.5;
decision = confHeads > 0.5;
ties = confHeads == 0.5;
decision(ties) = binornd(1,0.5,sum(ties),1);

%% P2: Psychometric conditioned on real block tendency
% For simplified plots of patterns (fitted)
x0 = [3 0.5 0.1 0.9];
lb = [0.1 -1.2 0 0.5];
ub = [30 1 0.5 1];
opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
ff = @(w,X) genlogf(X,w);

fObj = @(w) -f_obj( ff(w,mEv(blockBias==0)), confHeads(blockBias==0) );
w_T = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
p2.f_T = @(q) ff(w_T,q);

fObj = @(w) -f_obj( ff(w,mEv(blockBias==1)), confHeads(blockBias==1) );
w_H = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
p2.f_H = @(q) ff(w_H,q);

% q = 0.2:0.01:0.8;
% hold on
% plot(q,p2.f_T(q),'r');
% plot(q,p2.f_H(q),'b');

%% P4: pbEv match
% Align with real block bias
block = 2*(blockBias-0.5);
response = block.*(confHeads-0.5) + 0.5;
E = block.*(mEv-0.5) + 0.5;
P = block.*(pbEv-0.5) + 0.5;

% Simplified plot
linMod = fitlm(P,response,'linear');
p4.fx = @(q) sum([ones(size(q)) q].*linMod.Coefficients.Estimate',2);

% hold on
% plot(p4.x,p4.y);
% q = 0.3:0.01:0.9;
% plot(q',p4.fx(q'),'k');
% hold off

%% P6: sample crossover
% Align with real block bias
block = 2*(blockBias-0.5);
response = block.*(confHeads-0.5) + 0.5;
E = block.*(mEv-0.5) + 0.5;
P = block.*(pbEv-0.5) + 0.5;

many = [N>=quantile(N,0.6) N<=quantile(N,0.4)];

% For simplified plots of patterns (fitted)
x0 = [3 0.5 0.1 0.9];
lb = [0.1 -1.2 0 0.5];
ub = [30 1 0.5 1];
opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
ff = @(w,X) genlogf(X,w);

fObj = @(w) -f_obj( ff(w,E(many(:,1))), response(many(:,1)) );
w_H = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
p6.f_H = @(q) ff(w_H,q);

fObj = @(w) -f_obj( ff(w,E(many(:,2))), response(many(:,2)) );
w_T = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
p6.f_L = @(q) ff(w_T,q);

% hold on
% plot(p6.x,p6.y);
% q = 0.2:0.01:0.9;
% plot(q,p6.f_L(q),'g');
% plot(q,p6.f_H(q),'m');
% hold off

%% P7: prev. sample crossover
% Block evidence prior to present trial
pmEv = [nan; mEv(1:end-1)];
pmEv(1:blockLength:numel(mEv)) = nan;
pN = [nan; N(1:end-1)];
pN(1:blockLength:numel(mEv)) = nan;

% Align with real block bias
block = 2*(blockBias-0.5);
response = block.*(confHeads-0.5) + 0.5;
E = block.*(pmEv-0.5) + 0.5;

% Remove first in-block trial
idx = 1:blockLength:numel(mEv);
response(idx) = [];
E(idx) = [];
pN(idx) = [];

many = [pN>=quantile(pN,0.5) pN<quantile(pN,0.5)];
subset = any(many,2);

% For simplified plots of patterns (fitted)
x0 = [3 0.5 0.1 0.9];
lb = [0.1 -1.2 0 0.5];
ub = [30 1 0.5 1];
opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
ff = @(w,X) genlogf(X,w);

fObj = @(w) -f_obj( ff(w,E(many(:,1))), response(many(:,1)) );
w_H = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
p7.f_H = @(q) ff(w_H,q);

fObj = @(w) -f_obj( ff(w,E(many(:,2))), response(many(:,2)) );
w_T = fmincon(fObj,x0,[],[],[],[],lb,ub,[],opts);
p7.f_L = @(q) ff(w_T,q);

% hold on
% plot(p7.x,p7.y);
% q = 0.2:0.01:0.9;
% plot(q,p7.f_L(q),'g');
% plot(q,p7.f_H(q),'m');
% hold off

%% P8: Recency
% Bias seems to stem from different number of predictors
% The more the higher the weight
E = (mEv-0.5)*2;
% Dependent variable
mask = any(inBlock==[2:5],2);
expY = confHeads(mask);

% Index corresponds to previous trial
Xdat = zeros(numSim,blockLength-1);
numP = inBlock-1;
for j=1:4
    x = E(any(inBlock==[1:j],2));
    Xdat(any(inBlock==[5-(j-1):5],2),j) = x;
end
Xdat(inBlock==1,:) = [];
numP(inBlock==1,:) = [];

wExp = fit_prior_weights( fliplr(Xdat), expY, numP );
p8.x = -4:-1;
p8.y = fliplr(wExp);

%% P9: Accumulation of evidence across trials
% Align with real block bias
block = 2*(blockBias-0.5);
response = block.*(confHeads-0.5) + 0.5;
E = block.*(mEv-0.5) + 0.5;

bMask = inBlock==1:blockLength;
R = nan(numSim/blockLength,blockLength);
for j=1:blockLength
    R(:,j) = response(bMask(:,j));
end

p9.x = 1:blockLength;
p9.y = mean(R);

%% Plot
nXpanels = 3;
nYpanels = 2;
x_lmargin = 0.11;
x_rmargin = 0.05;
x_relWidth = [0.333 0.333 0.333];
x_spacing = 0.16;
y_lmargin = 0.12;
y_umargin = 0.1;
y_relWidth = [0.5 0.5];
y_spacing = 0.17;

% x and y extent
xe = (1-(x_lmargin+x_rmargin+2*x_spacing))*x_relWidth;
ye = (1-(y_lmargin+y_umargin+1*y_spacing))*y_relWidth; 
% Left edge
sep = [x_lmargin xe(1) x_spacing xe(2) x_spacing xe(3) x_rmargin];  % all separations
% assert(sum(sep)==1);
sep = cumsum(sep);
xedge = [sep(1) sep(3) sep(5)];
% Lower edge
sep = [y_lmargin ye(1) y_spacing ye(2) y_umargin];  % all separations
% assert(abs(sum(sep)-1)<=0.01);
sep = cumsum(sep);
yedge = [sep(1) sep(3)];

width = 17;
height = 11.33;
LW = 1.2;
FS = 9;

figure(1)
clf;


set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% P2 Psychometric
subplot('Position',[xedge(1) yedge(2) xe(1) ye(2)]);
hold on
line([0 100], [0.5 0.5], 'Color', [1 1 1]*0.8, 'LineWidth', 0.8);
line([50 50], [0 1], 'Color', [1 1 1]*0.8, 'LineWidth', 0.8);
q = [0.2:0.01:0.8]';
clear h
h(2) = plot(q*100, p2.f_T(q), 'Color', hsv2rgb([.86 1 0.9]), 'LineWidth', LW+2);
h(1) = plot(q*100, p2.f_H(q), 'Color', hsv2rgb([.5 1 0.9]), 'LineWidth', LW+2);
xlim([0 100]);
ylim([0 1]);
% Custom legend
xl = xlim;
yl = ylim;
xpos = 7;
ypos = 0.86;
rectangle('Position',[xpos ypos diff(xl)/8 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.5 1 0.9]));
ypos = 0.75;
rectangle('Position',[xpos ypos diff(xl)/8 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.86 1 0.9]));

xlabel('blue samples (\%)', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel('confidence blue majority', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top', 'XTick', 0:25:100, 'YTick', 0:0.25:1);


% P4 Match of inferred tendency
subplot('Position',[xedge(2) yedge(2) xe(2) ye(2)]);
hold on
line([0.5 0.5], [0 1], 'LineStyle', '--', 'Color', [1 1 1]*0.8, 'LineWidth', 0.8);
q = [0.2:0.01:0.9]';
plot(q, p4.fx(q), 'Color', hsv2rgb([0 .7 0.9]), 'LineWidth', LW+2);
xlim([0.2 1]);
ylim([0.45 0.9]);
xlabel('aligned inferred tendency', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel('aligned confidence', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top', 'XTick', 0.25:0.25:1, 'YTick', 0.5:0.1:0.9);


% P6 sample crossover
subplot('Position',[xedge(3) yedge(2) xe(3) ye(2)]);
hold on
line([0.5 0.5]*100, [0 1], 'Color', [1 1 1]*0.8, 'LineWidth', 0.8);
q = [0.2:0.01:0.9]';
h(2) = plot(q*100, p6.f_L(q), 'Color', hsv2rgb([.08 1 0.85]), 'LineWidth', LW+2);
h(1) = plot(q*100, p6.f_H(q), 'Color', hsv2rgb([.33 1 .85]), 'LineWidth', LW+2);
xlim([10 100]);
ylim([0.2 1]);
% Custom legend
xl = xlim;
yl = ylim;
xpos = 54;
ypos = 0.35;
rectangle('Position',[xpos ypos diff(xl)/8 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.33 1 .85]));
ypos = 0.25;
rectangle('Position',[xpos ypos diff(xl)/8 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.08 1 0.85]));

xlabel('aligned sample proportion (\%) ', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel('aligned confidence', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top', 'XTick', 0:25:100, 'YTick', 0.25:0.25:1);


% P7 previous sample crossover
subplot('Position',[xedge(1) yedge(1) xe(1) ye(1)]);
hold on
line([0.5 0.5]*100, [0 1], 'Color', [1 1 1]*0.8, 'LineWidth', 0.8);
q = [0.2:0.01:0.9]';
h(2) = plot(q*100, p7.f_L(q), 'Color', hsv2rgb([.08 1 0.85]), 'LineWidth', LW+2);
h(1) = plot(q*100, p7.f_H(q), 'Color', hsv2rgb([.33 1 .85]), 'LineWidth', LW+2);
xlim([10 100]);
ylim([0.5 0.85]);
% Custom legend
xl = xlim;
yl = ylim;
xpos = 54;
ypos = 0.565;
rectangle('Position',[xpos ypos diff(xl)/8 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.33 1 .85]));
ypos = 0.52;
rectangle('Position',[xpos ypos diff(xl)/8 diff(yl)/30],'EdgeColor','none','FaceColor',hsv2rgb([.08 1 0.85]));

xlabel('aligned proportion of previous sample (\%)', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel('aligned confidence', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top', 'XTick', 0:25:100, 'YTick', 0.5:0.1:0.8);


% P8 Recency
subplot('Position',[xedge(2) yedge(1) xe(2) ye(1)]);
hold on
h(1) = line([1 2], [1 1]*nanmean(p8.y(2,:)), 'Color', [0.5 0.5 0.5], 'LineWidth', LW+2);
h(2) = line([1 3], [1 1]*nanmean(p8.y(3,:)), 'Color', [0.5 0.5 0.5], 'LineWidth', LW+2);
h(3) = line([1 4], [1 1]*nanmean(p8.y(4,:)), 'Color', [0.5 0.5 0.5], 'LineWidth', LW+2);
scatter([1:2 1:3 1:4], [ones(1,2)*nanmean(p8.y(2,:)) ones(1,3)*nanmean(p8.y(3,:)) ones(1,4)*nanmean(p8.y(4,:))], ...
    'Marker', 'x', 'LineWidth', LW+1, 'MarkerEdgeColor', 'k');
% lh = legend(h,'T3','T4','T5');
% set(lh,'Box','off','FontSize',FS-2,'Location','best');
xlim([0.5 4.5]);
ylim([1.1 1.4]);
xlabel('trial in block', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'tex');
ylabel('influence of previous trials', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top', 'XTick', 1:4);


% P9 Accumulation of evidence
subplot('Position',[xedge(3) yedge(1) xe(3) ye(1)]);
hold on
plot(p9.x, p9.y, 'Color', hsv2rgb([0 .7 0.9]), 'LineWidth', LW+2);
xlim([0.5 5.5]);
ylim([0.6 0.75]);
xlabel('trial in block', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel('aligned confidence', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top', 'XTick', 1:5, 'YTick', 0.6:0.05:0.75);

%% Print
print(gcf, '-dpng', '-r400', '.\..\..\plots\exp2\pattern_overview.png');