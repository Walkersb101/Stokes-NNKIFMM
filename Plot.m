data = readmatrix("data.csv");

F(57) = struct('cdata',[],'colormap',[]);
for j=1:57
    clf;
    ax = axes();
    view(ax, 3);
    xlim(ax, [-3 3]);
    ylim(ax, [-3 3]);
    zlim(ax, [-3 3]);
    hold on;
    for i = 0:(size(data(2:end,:),1)/9)-1
        x0 = data((i*9)+2:(i*9)+4,j)';
        b1 = data((i*9)+5:(i*9)+7,j);
        b2 = data((i*9)+8:(i*9)+10,j);
        b3 = cross(b1,b2);
        B = [b1 b2 b3];
        points = prolateSpheroid(0.1,0.5,0.05);
        points = (B * points')';
        points = points+x0.*ones(size(points,1),1);
        plot3(points(:,1),points(:,2),points(:,3),'.')
    end
    hold off;
    F(j) = getframe(gcf);
end
fig = figure;
movie(fig,F,2)