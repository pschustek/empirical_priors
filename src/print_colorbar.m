function [] = print_colorbar(cmap)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = size(cmap,1);

xmap = interp1(linspace(0,1,N)',cmap,linspace(0,1,3*N)');

N = size(xmap,1);

pvec = [repmat(0,N,1) linspace(1-1/N,0,N)' repmat(0.2,N,1) repmat(1/N,N,1)];
pvec = [repmat(0,N,1) linspace(0,1-1/N,N)' repmat(0.2,N,1) repmat(1/N,N,1)];

for j=1:N
    rectangle('Position',pvec(j,:),'EdgeColor','none','FaceColor',xmap(j,:));
end

rectangle('Position',[0 0 0.2 1],'EdgeColor','k','LineWidth',0.5);

end

