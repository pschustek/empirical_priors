% Fig. 1a
clear all
close all

%% Prepare
addpath('.\..\..\src\');

rng(2350);

width = 7;
height = 0.6*width;

LW = 1.5;
FS = 11;

%% Plot
figure(1)
clf;

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top', 'DataAspectRatio', [1 1 1]);

hold on

xD = 3;
yD = 2;
offScaleX = 0.8;
offScaleY = 0.4;

rectangle('Position',[offScaleX*xD offScaleY*yD xD yD],'FaceColor','none','EdgeColor','k','LineWidth', 2);
rectangle('Position',[0 0 xD yD],'FaceColor','w','EdgeColor','k','LineWidth', 2);

% Sample
[gridX, gridY, color] = get_sample_positions(4, 6, rng());
draw_sample(xD/2, yD/2, gridX, gridY, 1, color, 0);

% Response bar
cH = opt_inf.basic_confH(4,6,4,4);
draw_responses( offScaleX*xD+xD/2, offScaleY*yD+yD/2, cH, 1.4)

xlim([-0.1 5.5]);
ylim([-0.1 offScaleY*yD+yD+0.1]);
% xlabel('trial in block', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
% ylabel('confidence matching tendency', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');

set(gca,'visible','off', 'Position', [0 0 1 1]);

%% Print
print(gcf, '-dpng', '-r400', 'basic_task_schematic.png');