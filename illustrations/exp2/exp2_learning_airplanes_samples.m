% Fig. 3
clear all
close all

%% Prepare
addpath('.\..\..\src\');

rng(826002652)

% Get data
exp2_learning_example;

%% Plot airplanes
width = 1.4;
height = 0.3*width;
FS = 9;
LW = 0.7;
xBounds = [-0.5-0.0 0.5+0.0];
yBounds = [-0.15-0.0 0.15+0.0];

for j=1:4

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
        'DataAspectRatio', [1 1 1], 'Visible','off');
    
    hold on
    draw_airplane( 0, 0, mu(j), 1);
    
    xlim(xBounds);
    ylim(yBounds);
    
    % Print
    print(gcf, '-dpng', '-r400', ['prior_airplane_' num2str(j) '.png']);
end


%% Plot samples
width = 1.3;
height = width;
FS = 9;
LW = 0.7;
xBounds = [-0.6-0.02 0.6+0.02];
yBounds = [-0.6-0.02 0.6+0.02];

for j=1:4

    figure(1);
    clf;
    
    hold on
    [gridX, gridY, color] = get_sample_positions(nH(j), N(j), rng());
    draw_sample(0, 0, gridX, gridY, 1, color, 1);
    
    xlim(xBounds);
    ylim(yBounds);
    
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
        'DataAspectRatio', [1 1 1], 'Visible','off');
    
    % Print
    print(gcf, '-dpng', '-r400', ['prior_sample_' num2str(j) '.png']);
end