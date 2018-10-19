function [ ax ] = plot_row( D )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

FS = 9;
LW = 1.2;
nSub = size(D,1);

xBounds = [0.4 size(D,2)+0.6];
yBounds = [-30 30];

% Bounds according to Hartley interpretation
dHart5X = [xBounds xBounds(2:-1:1)];   % "substantial"
dHart5Y = [1 1 -1 -1]*5;
dHart20X = [xBounds xBounds(2:-1:1)];   % "decisive"
dHart20Y = [1 1 -1 -1]*20;
dBGX = [xBounds xBounds(2:-1:1)];
dBGY = [[1 1]*yBounds(2) [1 1]*yBounds(1)];

hold on
fill(dBGX,dBGY,[0.9 1 0.9],'EdgeColor','none');
fill(dHart20X,dHart20Y,[1 1 0.9],'EdgeColor','none');
fill(dHart5X,dHart5Y,[1 0.9 0.9],'EdgeColor','none');

%set(gca,'XDir','reverse');
for j=1:size(D,2)
    % Individual differences
    Y = D(:,j);
    X = ones(1,nSub)*j;
    scatter(X,Y,'Marker','.','MarkerEdgeColor','k','SizeData',40);
    % Across subject trend
    scatter(j,median(Y),'r','SizeData',50,'LineWidth',1.5);
    
    % Number of significant participants
    ngeq = sum(Y>=20);
    nleq = sum(Y<=-20);
    if ngeq > 0
       text(j-0.35,25,num2str(ngeq),'FontSize',8,'FontName','Cambria',...
           'HorizontalAlignment','center','VerticalAlignment','middle','Color','k'); 
    end
    if nleq > 0
       text(j-0.35,-25,num2str(nleq),'FontSize',8,'FontName','Times',...
           'HorizontalAlignment','right','VerticalAlignment','middle','Color','k'); 
    end
end

xlim(xBounds);
ylim(yBounds);

%ylabel('$\Delta$ CVLL','Interpreter','latex')

ax = gca;
set(ax, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out',...
    'XMinorTick', 'on', 'YMinorTick', 'off', 'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top',...
    'YTick',[-20 -10 0 10 20],'XTick',1:size(D,2),'XTickLabel',{},...
    'YTickLabel',{'-20','','0','','20'});

end

