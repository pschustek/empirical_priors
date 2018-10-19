function [] = draw_responses( posX, posY, cursor, scale)
% centimeter units

% Inputs
xRange = 1;
yRange = 0.2;
dotSizeRef = 0.1;   % in cm (unscaled)
dotSize = 72*(dotSizeRef*scale/2.54);   % in pts = 1/72 inch

%% Settings
LW = 2;

%% Plot
% Outer frame
rectangle('Position',[[posX posY]-[xRange yRange]*scale/2 [xRange yRange]*scale],'LineWidth',LW);
% Center line
line([1 1]*posX,[-0.18 0.18]*scale+posY,'LineWidth',LW,'Color','k');
% Vertical response cursor
line([1 1]*(cursor-0.5)*xRange*scale+posX,[-0.15 0.15]*scale+posY,'LineWidth',LW,'Color',[1 0.6 .2]);
% Category indicator
scatter(posX+xRange*scale/2,-0.22*scale+posY,'MarkerFaceColor','b','MarkerEdgeColor','none','SizeData',dotSize^2);
scatter(posX-xRange*scale/2,-0.22*scale+posY,'MarkerFaceColor','r','MarkerEdgeColor','none','SizeData',dotSize^2);

end