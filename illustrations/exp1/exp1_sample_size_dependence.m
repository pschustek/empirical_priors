% Fig. 1b
clear all
close all

%% Make content
rng(23);
h0 = 4;
t0 = 4;

areaColT = hsv2rgb([1 1 .9]);
areaColH = hsv2rgb([.61 1 .9]);

% [HgridX,HgridY] = ndgrid([0.97:-0.05:0.87],[3.3:-0.5:2.3]);
[HgridY,HgridX] = ndgrid([0.5 0],linspace(0,1,5));
HgridX = HgridX(:);
HgridY = HgridY(:);

% [TgridX,TgridY] = ndgrid([0.97:-0.05:0.87],[1.8:-0.5:0.8]);
[TgridY,TgridX] = ndgrid(1,linspace(0,1,5));
TgridX = TgridX(:);
TgridY = TgridY(:);

%% 
width = 6.5;
height = 3.6;
width = 8;
height = width*0.5538;
LW = 1.2;
FS = 11;
dotSize = 3^2;

clf;
set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'OuterPosition', [4 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0.01 0.01 width-0.01 height-0.01]);

q = 0:0.01:1;
ybound = [0 3.7];

% % Height
% lb = [0.7183 0.4187 0.1191];
% ve = 0.2067;
% % Width
% xb = [0.6916 0.4108 0.13];
% xe = 0.2134;

npanels = 3;
x_lmargin = 0.03;
x_rmargin = 0.03;
x_relWidth = [0.15 0.62 0.23];
x_spacing = 0.09;
y_lmargin = 0.16;
y_umargin = 0.03;
y_relWidth = [0.333 0.333 0.333];
y_spacing = 0.07;

% x and y extent
xe = (1-(x_lmargin+x_rmargin+2*x_spacing))*x_relWidth;
ye = (1-(y_lmargin+y_umargin+2*y_spacing))*y_relWidth; 
% Left edge
sep = [x_lmargin xe(1) x_spacing xe(2) x_spacing xe(3) x_rmargin];  % all separations
assert(sum(sep)==1);
sep = cumsum(sep);
xedge = [sep(1) sep(3) sep(5)];
% Lower edge
sep = [y_lmargin ye(1) y_spacing ye(2) y_spacing ye(3) y_umargin];  % all separations
assert(abs(sum(sep)-1)<=0.01);
sep = cumsum(sep);
yedge = [sep(1) sep(3) sep(5)];

%% First sample
H = 2;
T = 1;
% Sample schematic
%subplot(3,3,1);
%subplot('Position',[0.05 lb(1) 0.15 ve]);
subplot('Position',[xedge(1) yedge(3) xe(1) ye(3)]);
hold on
scatter(HgridX(1:H),HgridY(1:H),'MarkerFaceColor',areaColH,'SizeData',dotSize);
scatter(TgridX(1:T),TgridY(1:T),'MarkerFaceColor',areaColT,'SizeData',dotSize);
xlim([0 1]);
ylim([-0.25 1.25]);
axis off

% Distribution schematic
%subplot(3,3,2);
%subplot('Position',[0.23 lb(1) 0.45 ve]);
subplot('Position',[xedge(2) yedge(3) xe(2) ye(3)]);
hold on
xGrid = linspace(0,0.5,100);
yGrid = betapdf(xGrid,h0+H,t0+T);
xArea = [xGrid(1) xGrid xGrid(end)];
yArea = [0 yGrid 0];
fill(xArea,yArea,areaColT,'EdgeColor','none');
xGrid = linspace(0.5,1,100);
yGrid = betapdf(xGrid,h0+H,t0+T);
xArea = [xGrid(1) xGrid xGrid(end)];
yArea = [0 yGrid 0];
fill(xArea,yArea,areaColH,'EdgeColor','none');
plot(q, betapdf(q,h0+H,t0+T), 'k', 'LineWidth', LW);
ylim(ybound);
%xlabel('estimate');
set(gca,'YTick',[],'YTickLabel',{},'Layer', 'top', 'XTickLabel',{}, 'FontSize', FS);

% Histogram schematic
%subplot(3,3,3)
%subplot('Position',[0.8 lb(1) 0.15 ve]);
subplot('Position',[xedge(3) yedge(3) xe(3) ye(3)]);
hold on
cH = integral(@(q)betapdf(q,h0+H,t0+T),0.5,1);
hold on
bar(1,1-cH,'FaceColor',areaColT,'EdgeColor','none');
bar(2,cH,'FaceColor',areaColH,'EdgeColor','none');
xlim([0.5 2.5])
%xlabel('confidence');
set(gca,'XTick',[1 2],'XTickLabel',{}, 'YTick', [0 1],'YTickLabel',{}, 'FontSize', FS);

%% Second sample
H = 4;
T = 2;
% Sample schematic
%subplot(3,3,4);
%subplot('Position',[0.05 lb(2) 0.15 ve])
subplot('Position',[xedge(1) yedge(2) xe(1) ye(2)]);
hold on
scatter(HgridX(1:H),HgridY(1:H),'MarkerFaceColor',areaColH,'SizeData',dotSize);
scatter(TgridX(1:T),TgridY(1:T),'MarkerFaceColor',areaColT,'SizeData',dotSize);
xlim([0 1]);
ylim([-0.25 1.25]);
axis off

% Distribution schematic
% subplot(3,3,5);
% subplot('Position',[0.23 lb(2) 0.45 ve]);
subplot('Position',[xedge(2) yedge(2) xe(2) ye(2)]);
hold on
xGrid = linspace(0,0.5,100);
yGrid = betapdf(xGrid,h0+H,t0+T);
xArea = [xGrid(1) xGrid xGrid(end)];
yArea = [0 yGrid 0];
fill(xArea,yArea,areaColT,'EdgeColor','none');
xGrid = linspace(0.5,1,100);
yGrid = betapdf(xGrid,h0+H,t0+T);
xArea = [xGrid(1) xGrid xGrid(end)];
yArea = [0 yGrid 0];
fill(xArea,yArea,areaColH,'EdgeColor','none');
plot(q, betapdf(q,h0+H,t0+T), 'k', 'LineWidth', LW);
ylim(ybound);
ylabel('pdf');
set(gca,'YTick',[],'YTickLabel',{},'Layer', 'top', 'XTickLabel',{}, 'FontSize', FS);

% Histogram schematic
% subplot(3,3,6)
% subplot('Position',[0.8 lb(2) 0.15 ve]);
subplot('Position',[xedge(3) yedge(2) xe(3) ye(2)]);
hold on
cH = integral(@(q)betapdf(q,h0+H,t0+T),0.5,1);
hold on
bar(1,1-cH,'FaceColor',areaColT,'EdgeColor','none');
bar(2,cH,'FaceColor',areaColH,'EdgeColor','none');
xlim([0.5 2.5])
ylabel('probability');
set(gca,'XTick',[1 2],'XTickLabel',{}, 'YTick', [0 1],'YTickLabel',{}, 'FontSize', FS);

%% Third sample
H = 8;
T = 4;

% Sample schematic
% subplot(3,3,7);
% subplot('Position',[0.05 lb(3) 0.15 ve]);
subplot('Position',[xedge(1) yedge(1) xe(1) ye(1)]);
hold on
scatter(HgridX(1:H),HgridY(1:H),'MarkerFaceColor',areaColH,'SizeData',dotSize);
scatter(TgridX(1:T),TgridY(1:T),'MarkerFaceColor',areaColT,'SizeData',dotSize);
xlim([0 1]);
ylim([-0.25 1.25]);
lh = xlabel('sample');
lh.Units = 'data';
pos = lh.Position(1:2);
xlabel('sample', 'FontSize', FS);
axis off
text(gca,pos(1),pos(2),'sample','HorizontalAlignment','center','VerticalAlignment','cap', 'FontSize', FS, 'Units', 'data');

% Distribution schematic
% subplot(3,3,8);
% subplot('Position',[0.23 lb(3) 0.45 ve]);
subplot('Position',[xedge(2) yedge(1) xe(2) ye(1)]);
hold on
xGrid = linspace(0,0.5,100);
yGrid = betapdf(xGrid,h0+H,t0+T);
xArea = [xGrid(1) xGrid xGrid(end)];
yArea = [0 yGrid 0];
fill(xArea,yArea,areaColT,'EdgeColor','none');
xGrid = linspace(0.5,1,100);
yGrid = betapdf(xGrid,h0+H,t0+T);
xArea = [xGrid(1) xGrid xGrid(end)];
yArea = [0 yGrid 0];
fill(xArea,yArea,areaColH,'EdgeColor','none');
plot(q, betapdf(q,h0+H,t0+T), 'k', 'LineWidth', LW);
ylim(ybound);
% xlabel('estimate', 'FontSize', FS);
set(gca,'YTick',[],'YTickLabel',{},'Layer', 'top', 'FontSize', FS, 'XTick', 0:0.5:1,'XTickLabel',{}, 'TickDir','out');

% Histogram schematic
% subplot(3,3,9)
% subplot('Position',[0.8 lb(3) 0.15 ve]);
subplot('Position',[xedge(3) yedge(1) xe(3) ye(1)]);
hold on
cH = integral(@(q)betapdf(q,h0+H,t0+T),0.5,1);
hold on
bar(1,1-cH,'FaceColor',areaColT,'EdgeColor','none');
bar(2,cH,'FaceColor',areaColH,'EdgeColor','none');
xlim([0.5 2.5])
% hl = xlabel('confidence', 'FontSize', FS);
set(gca,'XTick',[1 2],'XTickLabel',{}, 'YTick', [0 1],'YTickLabel',{}, 'FontSize', FS);

%% Print
print(gcf, '-dpng', '-r400', 'sample_size_dependence.png');