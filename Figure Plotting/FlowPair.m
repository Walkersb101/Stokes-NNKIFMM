% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder))

Evalres = 50;
[X,Y,Z] = meshgrid(linspace(-2,2,Evalres),[-1 0 1],linspace(-1,3,Evalres));
EvalPoints = [X(:) Y(:) Z(:)];

raw = readmatrix('SquiremerNearPair.csv');

coarse = prolateSpheroid(0.5,0.5,0.2);
velsquirmer = squirmerVel(coarse,1.5);
kernelPar = [1e-2,1];

fine = prolateSpheroid(0.5,0.5,0.05);
fine = unique([fine;coarse],'rows','stable');
NN = nearestMatrix(coarse,fine,1,0);


framesSteps = linspace(min(raw(1,:)),max(raw(1,:)),8);
data = spline(raw(1,:),raw(2:end,:),framesSteps);

noFrames = size(data,2);
for j=1:noFrames
    figure('Position', [0 0 1100 1000])
    ax = axes();
    view(ax, 3);
    view(0,0)
    hold on;
    
    squirmer = Swimmers();
    initialCon = [];
    
%     bndCoarse = [wallTri([0 0 0],[10 10],3,0.5);...
%              wallTri([0 0 10],[10 10],3,0.5)];
%     squirmer.updateBnd(bndCoarse,zeros(size(bndCoarse)));

    
    for i = 0:(size(data,1)/9)-1
        x0 = data((i*9)+1:(i*9)+3,j)';
        b1 = data((i*9)+4:(i*9)+6,j);
        b2 = data((i*9)+7:(i*9)+9,j);
        
        squirmer.addSwimmer(coarse,fine,NN,x0,[b1;b2],...
                    velsquirmer,[0;0;0],[0;0;0]);
                
        initialCon = vertcat(initialCon,x0',[b1;b2]);
    end
    
    [points,fine2,NN2] = squirmer.getSwimmerBnd();
    
    FMMEval = KIFMM(OcTree(fine2,EvalPoints),kernelPar,'parThreads',20,'GMRES',1,'format',2);
    
    shearFunc = @(t,x) shearFlow(x,0);
    kernelPar = [1e-2,1];
    treePar = {'nodeCapacity', 500, 'nearest', 1};
    
    pot = mobilityProblemGravityForce(squirmer, shearFunc,0,initialCon,...
                                      kernelPar,treePar,20,0,1,1e-4,75,...
                                      [2e-5,100],[0;0;0;0]);
                                  
    pot = kron(eye(3),NN2)*pot;
    
    EvalVel = reshape(FMMEval.computeVel(pot),[],3);
    vx = reshape(EvalVel(:,1),[],Evalres,Evalres);
    vy = reshape(EvalVel(:,2),[],Evalres,Evalres);
    vz = reshape(EvalVel(:,3),[],Evalres,Evalres);
    
    velMag = sqrt(vx.^2 + vy.^2 + vz.^2);
    s = slice(X,Y,Z,velMag,[],0,[]);
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    
    for i = 0:(size(data,1)/9)-1
        x0 = data((i*9)+1:(i*9)+3,j)';
        b1 = data((i*9)+4:(i*9)+6,j);
        b2 = data((i*9)+7:(i*9)+9,j);
        b3 = cross(b1,b2);
        B = [b1 b2 b3];
        points = prolateSpheroid(0.5,0.5,0.2);
        points = (B * points')';
        points = points+x0.*ones(size(points,1),1);
        plot3(ax,points(:,1),points(:,2),points(:,3),'w.','MarkerSize',20)
        points = prolateSpheroid(0.5,0.5,0.01);
        points = (B * points')';
        points = points+x0.*ones(size(points,1),1);
        plot3(ax,points([1 end],1),points([1 end],2),points([1 end],3),'w','LineWidth',6)
    end
    
%     [startX,startY,startZ] = meshgrid(-2.5:0.25:2.5,[0],[-1 4]);
%     verts = stream3(X,Y,Z,vx,vy,vz,startX,startY,startZ,[0.005,50000]);
%     l = streamline(verts);
%     set(l,'Color','w');
    
    q = quiver3(X(3,1:4:end,1:4:end),zeros(size(X(3,1:4:end,1:4:end))),Z(3,1:4:end,1:4:end),vx(3,1:4:end,1:4:end),zeros(size(X(3,1:4:end,1:4:end))),vz(3,1:4:end,1:4:end));
%     q.AutoScaleFactor = 2;
    q.Color = 'w';
%     q.LineWidth = 1;
%     q.MarkerSize = 10;
%     q.MaxHeadSize = 20;
    
    zlabel(ax,"x_3", 'FontSize',40);
    xlabel(ax,"x_2", 'FontSize', 40);
    h = colorbar('FontSize',40);
    set(gca,"FontSize",40);
    set(gca,'ColorScale','log');
    ylabel(h, 'Magnitude of velocity', 'FontSize', 40)
%     caxis([0 1.5])
    
    xlim(ax, [-2 2]);
    ylim(ax, [-2 2]);
    zlim(ax, [-1 3]);
    hold off;
    
    daspect(ax,[1 1 1])
    ax.Position = [0.1529    0.1507    0.62    0.7743];
    h.Position = [0.778    0.196    0.0500    0.683];
    
    savefig(gcf,sprintf('Pair-%d.pdf',j))
end
display(frameSteps)

function vel = squirmerVel(points,scale)
tangent = [points(:,1).*points(:,3), points(:,2).*points(:,3),...
           -(points(:,1).^2 + points(:,2).^2)]...
           ./(sqrt(points(:,1).^2 + points(:,2).^2+ points(:,3).^2)...
           .*sqrt(points(:,1).^2 + points(:,2).^2));
tangent(isnan(tangent))=0;
theta = acos(points(:,3)./...
        (sqrt(points(:,1).^2 + points(:,2).^2 + points(:,3).^2)));
    
vel = scale*sin(theta).*tangent;
end

function vel = shearFlow(points, gradient)
vel = zeros(size(points));
vel(:,1) = gradient*points(:,3);
end

function savefig(fig,name)
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,name,'-dpdf','-r0')
end