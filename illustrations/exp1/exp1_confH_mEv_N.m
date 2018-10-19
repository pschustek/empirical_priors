% Fig. 1c
clear all
close all

%% Prepare
addpath('.\..\..\src\')

h0 = 4;
t0 = 4;

%% Plot
width = 8;
height = 5.5;
FS = 11;
LW = 1.2;
figure(1)
clf;
hold on

colGrad = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*1, ones(numel(q),1)*0.85]);

vec_N = 3:13;

line([1 1]*50,[0 1],'LineStyle','--','LineWidth',0.8,'Color',[1 1 1]*0.8);
line([0 100],[1 1]*0.5,'LineStyle','--','LineWidth',0.8,'Color',[1 1 1]*0.8);

j = 0;
for N=vec_N
    j = j + 1;
    legLabel{j} = ['N = ' num2str(vec_N(j))];

    H = 0:N;
    T = N - H;
    confH = nan(numel(N),1);
    for n=1:numel(H)
        confH(n) = integral(@(q)betapdf(q,h0+H(n),t0+T(n)),0.5,1);
    end
    h(j) = plot(H/N*100,confH,'Color',colGrad((j-1)/(numel(vec_N)-1)),'LineWidth',LW);
end

xlim([0 100]);
%ylim([0.5 1]);
xlabel('proportion of blue samples (\%)', 'FontSize', FS, 'FontName', 'Times','Interpreter','latex');
ylabel('confidence blue majority', 'FontSize', FS, 'FontName', 'Times','Interpreter','latex');

% Color gradients with patch()
% lh = legend(h([1 end]), legLabel{[1 end]});
% set(lh, 'Location', 'northwest', 'FontSize', FS-1, 'FontName', 'Times','box','off');

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top','XTick',0:20:100);

%% Print
print(gcf, '-dpng', '-r400', 'confH_sample_size.png');

%% Plot colorbar separately
colGrad = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*1, ones(numel(q),1)*0.85]);
cmap = colGrad(0:0.02:1);

cbh = 3;
cbw = cbh*0.2;
clf;
set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 cbw cbh]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [cbw cbh], 'PaperPosition', [0 0 cbw cbh]);

print_colorbar(cmap);
set(gca,'DataAspectRatio',[1 1 1], 'Box', 'off', 'FontSize', 9, 'FontName', 'Times', 'TickDir', 'out', ...
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top', 'InnerPosition', [0 0 1 1]);
axis off
print(gcf, '-dpng', '-r400', 'sample_size_colorbar.png');