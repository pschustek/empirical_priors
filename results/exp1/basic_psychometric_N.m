%% Binned plot of N dependence of psychometric function
% Fig. 2a
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\src\');

% Subselect subjects
subInd = 1:24;

% Load data
load('.\..\..\data\exp1_data.mat');

%% Compute results
setN = [3 5 7 8 9 10 11 12 13]';

for n=1:length(setN)    
    % Mask
    H = 0:setN(n);
    T = setN(n) - H;
    
    bin{n}.x = H./setN(n);
    % N-specific bins
    for b=1:length(H)
        R = [];
        for s=subInd
            clear trials
            trials = trialData{s};
            
            mEv = trials.meanEvidence;
            N = trials.sampleSize;
            response = trials.confHeads;
            
            assert(mean(unique(N)==setN)==1);
    
            nMask = N == setN(n);
            mEvMask = mEv == H./(H+T);
            mask = nMask & mEvMask;
            
            R = [R; response(mask(:,b))];
        end
       
        bin{n}.mean(b,1) = mean(R);
        bin{n}.median(b,1) = median(R);
        bin{n}.std(b,1) = std(R);
        bin{n}.SE(b,1) = std(R)/sqrt(length(R));
    end
    
end


%% Optimal behavior
h0 = 4;
t0 = 4;
for n=1:length(setN)
    
    H = 0:setN(n);
    T = setN(n) - H;
    
    opt{n}.x = H./setN(n);
    opt{n}.conf = opt_inf.basic_confH( H, H+T, h0, t0 );
end

%% Plot
width = 8;
height = 6;
LW = 1.2;
FS = 11;
figure(1)
clf;
hold on

colGrad = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*1, ones(numel(q),1)*0.85]);
colGradModel = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*.35, ones(numel(q),1)*1]);


line([1 1]*50,[0 1],'LineStyle','--','LineWidth',0.8,'Color',[1 1 1]*0.8);
line([0 100],[1 1]*0.5,'LineStyle','--','LineWidth',0.8,'Color',[1 1 1]*0.8);

for n=1:length(setN)
    plot(opt{n}.x*100,opt{n}.conf,'Color',colGradModel((n-1)/(numel(setN)-1)),'LineWidth',LW+2,'LineStyle','-');
end

for n=1:length(setN)
    h(n) = errorbar(bin{n}.x*100,bin{n}.mean,bin{n}.SE,'Color',colGrad((n-1)/(numel(setN)-1)),'LineWidth',LW,'LineStyle','-','CapSize',0);
end

xlim([0 100]);
%ylim([0.5 1]);

xlabel('proportion of blue samples (\%)', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel('confidence blue majority', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');

%% Print
print(gcf, '-dpng', '-r400', '.\..\..\plots\exp1\basic_psychometric_N.png');
