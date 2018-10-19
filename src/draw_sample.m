function [] = draw_sample(posX, posY, gridX, gridY, scale, color, frame)
% for unit grid

% gridX = gridX - 0.5;
% gridY = gridY - 0.5;

% Inputs
% maxN = 13;
% dotSize = 2/(1.6*(maxN-1)), according to generation
dotSizeRef = .15;   % in cm (unscaled)
dotSize = 72*(dotSizeRef*scale/2.54);   % in pts = 1/72 inch
LW = 1.8*scale;
frameOffset1 = 0.2;   % offset as factor
frameOffset2 = 0.5;

range = scale*(1+frameOffset1);
frame1 = [posX-range/2 posY-range/2 range range];
range = scale*(1+frameOffset2);
frame2 = [posX-range/2 posY-range/2 range range];

switch frame
    case 1
        % Add one surrounding black frame
        rectangle('Position',frame1,'LineWidth',LW);
    case 2
        % Add one surrounding black frame
        rectangle('Position',frame2,'LineWidth',LW,'FaceColor',[0.8 0.1 0.1],'EdgeColor','none');
        % Add second surrounding frame 
        rectangle('Position',frame1,'LineWidth',LW,'FaceColor','w');
    case 3
        % Add one surrounding black frame
        rectangle('Position',frame2,'LineWidth',LW,'FaceColor',[0.1 0.8 0.1],'EdgeColor','none');
        % Add second surrounding frame
        rectangle('Position',frame1,'LineWidth',LW,'FaceColor','w');
    otherwise
end
% % Add one surrounding black frame
% rectangle('Position',[posX-scale*frame2/2 posY-scale*frame2/2 scale*(1+frame2) scale*(1+frame2)],'LineWidth',LW,'FaceColor',[0.8 0.1 0.1]);
% % Add second surrounding frame 
% rectangle('Position',[posX-scale*frame1/2 posY-scale*frame1/2 scale*(1+frame1) scale*(1+frame1)],'LineWidth',LW,'FaceColor','w');

%% Plot 
dotPos = [posX+scale*gridX; posY+scale*gridY];

scatter(dotPos(1,:),dotPos(2,:),dotSize.^2,color','filled')

end