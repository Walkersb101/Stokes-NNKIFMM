v = VideoWriter('video.avi');
open(v);

data = readmatrix("SquiremerGrav.csv");

squiremer = prolateSpheroid(0.5,0.5,0.02);

noFrames = size(data,2);
for j=1:noFrames
    set(gcf,'position',[0,0,1000,1000])
    clf;
    ax = axes();
    view(ax, 3);
    view(45,45)
    xlim(ax, [-5 5]);
    ylim(ax, [-5 5]);
    zlim(ax, [-5 5]);
    hold on;
    for i = 0:(size(data(2:end,:),1)/9)-1
        x0 = data((i*9)+2:(i*9)+4,j)';
        b1 = data((i*9)+5:(i*9)+7,j);
        b2 = data((i*9)+8:(i*9)+10,j);
        b3 = cross(b1,b2);
        B = [b1 b2 b3];
        points = squiremer;
        points = (B * points')';
        points = points+x0.*ones(size(points,1),1);
        plot3(points(:,1),points(:,2),points(:,3),'.')
    end
    hold off;
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);