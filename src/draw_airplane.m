function [] = draw_airplane( posX, posY, frac, scale)
% centimeter units

% Inputs
xRange = 1;
yRange = 0.3;

%% Settings
LW = 1.8*scale;

%% Plot
% Red 
rectangle('Position',[posX-scale*xRange/2 posY-scale*yRange/2 scale*xRange scale*yRange],'FaceColor','r','EdgeColor','none');
% Blue (partial)
rectangle('Position',[posX-scale*xRange/2+scale*xRange*(1-frac) posY-scale*yRange/2 scale*xRange*frac scale*yRange],'FaceColor','b','EdgeColor','none');

% Frame
rectangle('Position',[posX-scale*xRange/2 posY-scale*yRange/2 scale*xRange scale*yRange],'FaceColor','none','EdgeColor','k','LineWidth',LW);

end