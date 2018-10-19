% Fig. 3
clear all
close all

%% Prepare
addpath('.\..\..\src\');

% Parameters
width = 1.5;
height = 0.6*width;

h0 = 14;
t0 = 9;

colPos = [0 0.8 0.8];
colNeg = [0.8 0 0.8];

bgCol = [220 255 220]/255;

colPosArea = [0.7 1 1];
colNegArea = [1 0.7 1];

mode = (h0-1)/(h0+t0-2);
maxpdf = betapdf(mode,h0,t0);

xOffPdf = 1;
yOffPdf = 1;

scl = width/height;

%% Plot
FS = 9;
LW = 0.7;
xBounds = [0-0.05 1+0.05] + xOffPdf;
yBounds = [0.9 1/scl+yOffPdf+0.05];

figure(1)
hold on

x = linspace(0,1,500);

clf;
clear h
hold on

off = 0.95;
draw_airplane( 1, off, 0, 0.15);
draw_airplane( 1.5, off, .5, 0.15);
draw_airplane( 2, off, 1, 0.15);

% Area corresponding depending on skewness type
xGrid = linspace(0,1,200);      
yGrid = betapdf(xGrid,h0,t0)/(scl*maxpdf);
xArea = [xGrid(1) xGrid xGrid(end)] + xOffPdf;
yArea = [0 yGrid 0] + yOffPdf;
if h0 > t0
    fill(xArea,yArea,colPosArea, 'EdgeColor','none');
    h(1) = plot(x+xOffPdf,betapdf(x,h0,t0)/(scl*maxpdf) + yOffPdf,'Color',colPos,'LineWidth',LW+0.4);
else
    fill(xArea,yArea,colNegArea, 'EdgeColor','none');
    h(1) = plot(x+xOffPdf,betapdf(x,h0,t0)/(scl*maxpdf) + yOffPdf,'Color',colNeg,'LineWidth',LW+0.4);
end

line([1 1]*0.5 + xOffPdf,[0 betapdf(0.5,h0,t0)/(scl*maxpdf)] + yOffPdf,'Color','k','LineWidth',LW);
line([0 1]+xOffPdf,[1 1]*yOffPdf,'Color','k','LineWidth',LW);

 
xlim(xBounds);
ylim(yBounds);

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
    'XMinorTick', 'on', 'YMinorTick', 'off', 'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top', ...
    'XTick',[xBounds(1):0.5:xBounds(2)],'YTick',[],'YTickLabel',{}, 'DataAspectRatio', [1 1 1],...
    'Visible','off');

%% Print
print(gcf, '-dpng', '-r400', ['beta_' num2str(h0) '_' num2str(t0) '_w' num2str(round(width*10)) 'mm_' '.png']);