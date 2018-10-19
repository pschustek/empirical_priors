% Fig. 3
clear all
close all

%% Prepare
addpath('.\..\..\src\');

h0 = 14;
t0 = 9;
mode = (h0-1)/(h0+t0-2);
maxpdf = betapdf(mode,h0,t0);

% Get data
exp2_learning_example;

% Plot bottom-up message for first 3 trial and 
% top-down belief about block tendency M(b) for 4th trial
frac = [bEv(1:3); pbEv(4)];

%% Plot airplanes
width = 1.1;
height = 1.3*width;
xD = 0.8;   % bar width
labelOffset = -0.3;     % switch on skewed betas
bgCol = [255 240 220]/255; % [1 1 1]
FS = 9;
LW = 1.2;
xBounds = [-0.5-0.02 1.5+0.02];
yBounds = [-0.08+labelOffset 1+0.02];

for j=1:4

    figure(1);
    clf;
    
    % Plots
    set(gcf,'Color',bgCol,'InvertHardcopy','off'); 
    
    % Position plot on the screen for drawing
    set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);
    
    % Position plot on the paper for printing
    set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
        'PaperSize', [width height], 'PaperPosition', [0 0 width height]);
    
    % Axes
    % 'Position' positions drawing pane within window [0.15 0.15 0.78 0.8] for FS=18
    % 'OuterPosition' positions drawing pane together WITH axis within window
    set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'Position', [0 0 1 1], ...
        'XMinorTick', 'off', 'YMinorTick', 'off', 'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top', ...
        'XTick',1:2, 'YTick',0:0.5:1, 'XTickLabel',{}, 'YTickLabel',{});
    
    hold on
    rectangle('Position',[-xD/2 0 xD 1-frac(j)], 'EdgeColor','none', 'FaceColor','m','LineWidth',LW);
    rectangle('Position',[-xD/2+1 0 xD frac(j)], 'EdgeColor','none', 'FaceColor','c','LineWidth',LW);
    rectangle('Position',[-0.5 0 2 1], 'EdgeColor','k', 'FaceColor','none','LineWidth',LW);
    
    q = linspace(0,1,100);
    downsize = 1.3;
    shiftFac = 1.8;
    
    yBeta = betapdf(q,t0,h0)/maxpdf*0.2;
    plot(q/downsize-xD/shiftFac,yBeta+labelOffset,'Color','m','LineWidth',1);
    line([0 1]/downsize-xD/shiftFac,[1 1]*labelOffset,'Color','k');
    line([.5 .5]/downsize-xD/shiftFac,[1 .2]*labelOffset,'Color','k');
    
    yBeta = betapdf(q,h0,t0)/maxpdf*0.2;
    plot(q/downsize-xD/shiftFac+1,yBeta+labelOffset,'Color','c','LineWidth',1);
    line([0 1]/downsize-xD/shiftFac+1,[1 1]*labelOffset,'Color','k');
    line([.5 .5]/downsize-xD/shiftFac+1,[1 .2]*labelOffset,'Color','k');
    
    xlim(xBounds);
    ylim(yBounds);
    
    % Print
    print(gcf, '-dpng', '-r400', ['prior_messages_' num2str(j) '.png']);
end

%% Plot response on last trial
width = 2;
height = 0.4*width;
scl = 2;
xBounds = [-0.5 0.5]*scl + [-1 1]*0.1*scl;
yBounds = [-0.5/2 0.3/2]*scl + [-1 1]*0.05*scl;
figure(1);
clf;

% Plots
set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
% 'Position' positions drawing pane within window [0.15 0.15 0.78 0.8] for FS=18
% 'OuterPosition' positions drawing pane together WITH axis within window
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'Position', [0 0 1 1], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top', ...
    'Visible','off', 'DataAspectRatio', [1 1 1]);

hold on
draw_responses( 0, 0, confHeads(4), scl)

xlim(xBounds);
ylim(yBounds);

% Print
print(gcf, '-dpng', '-r400', ['response_trial_4.png']);
