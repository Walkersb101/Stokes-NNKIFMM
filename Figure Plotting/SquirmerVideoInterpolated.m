% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder))


v = VideoWriter('videoInter.avi');
open(v);

raw = readmatrix('SquiremerNearPair.csv');

squiremer = prolateSpheroid(0.5,0.5,0.05);

framesSteps = min(raw(1,:)):(1/60):max(raw(1,:));
data = spline(raw(1,:),raw(2:end,:),framesSteps);

noFrames = size(data,2);
for j=1:noFrames
    clf;
    set(gcf,'position',[100,100,400,400])
    ax = axes();
    view(ax, 3);
    view(0,0)
    hold on;
    for i = 0:(size(data,1)/9)-1
        x0 = data((i*9)+1:(i*9)+3,j)';
        b1 = data((i*9)+4:(i*9)+6,j);
        b2 = data((i*9)+7:(i*9)+9,j);
        b3 = cross(b1,b2);
        B = [b1 b2 b3];
        points = squiremer;
        points = (B * points')';
        points = points+x0.*ones(size(points,1),1);
        plot3(ax,points(:,1),points(:,2),points(:,3),'.')
    end
    xlim(ax, [-2 2]);
    ylim(ax, [-2 2]);
    zlim(ax, [-2 2]);
    hold off;
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);