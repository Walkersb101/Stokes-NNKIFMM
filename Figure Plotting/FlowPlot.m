[startX,startY,startZ] = meshgrid(-2:0.1:2,0,2);

evalVel = readmatrix("singleFlow.csv");
Evalres = 50;
[X,Y,Z] = meshgrid(linspace(-3,3,Evalres),linspace(-3,3,Evalres),linspace(-3,3,Evalres));


vx = reshape(evalVel(:,1),Evalres,Evalres,Evalres);
vy = reshape(evalVel(:,2),Evalres,Evalres,Evalres);
vz = reshape(evalVel(:,3),Evalres,Evalres,Evalres);

verts = stream3(X,Y,Z,vx,vy,vz,startX,startY,startZ);                
                  
velmag = sqrt(vx.^2 + vy.^2 + vz.^2);

figure('Position', [200 200 800 600])
hold on
% [startX,startY] = meshgrid(-2:0.25:2,2);
% stream2(x,y,u,v,startx,starty)
l = streamline(verts);
set(l,'Color','w');
view(0,0)
% plot3(coarse(:,1),coarse(:,2),coarse(:,3),'b.');
s = slice(X,Y,Z,velmag,[],0,[],'nearest');
s.FaceColor = 'interp';
s.EdgeColor = 'non';

zlabel(gca,"x_3",'FontSize',24);
xlabel(gca,"x_2",'FontSize',24);
h = colorbar('Location','eastoutside');
ylabel(h, 'Magnitude of velocity')
axis equal;
xlim([-1.5,1.5])
ylim([-1.5,1.5])
zlim([-1.5,1.5])
hold off
set(gca,"FontSize",24)
set(gca,'ColorScale','log')
% h.Position = [0.75    0.2509    0.0492    0.6741];
savefig(gcf,'StreamLinesSingle.pdf')


function savefig(fig,name)
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,name,'-dpdf','-r0')
end