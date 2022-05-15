clear all; close all; clc;

spmd
   gpuDevice(1);
end

squirmer = Swimmers();
initialCon = [];

bndCoarse = vertcat(wall([5 0 0],[10 10],1,0.3),...
                   wall([-5 0 0],[10 10],1,0.3),...
                   wall([0 5 0],[10 10],2,0.3),...
                   wall([0 -5 0],[10 10],2,0.3),...
                   wall([0 0 5],[10 10],3,0.3),...
                   wall([0 0 -5],[10 10],3,0.3));
               
bndCoarse = unique(bndCoarse,'rows','stable');
               
NN = nearestMatrix(bndCoarse,bndCoarse,1,0);
squirmer.updateBnd(bndCoarse,bndCoarse,NN,zeros(size(bndCoarse)));
bndFine = []; bndCoarse = [];

coarse = prolateSpheroid(0.5,0.5,0.25);
fine = prolateSpheroid(0.5,0.5,0.13);
NN = nearestMatrix(coarse,fine,1,0);
vel = squirmerVel(coarse,2);

startlocal = -4 + 8*rand(10,3);

initialCon = [];
for i = 1:10
    basis = Rotate(eye(3),2*pi*rand(1),2*pi*rand(1),2*pi*rand(1),[0 0 0]);
    b1b2 = basis(1:6)';
    
    squirmer.addSwimmer(coarse,fine,NN,startlocal(i,:),b1b2,...
                        vel,[0;0;0],[0;0;0]);
    initialCon = vertcat(initialCon,startlocal(i,:)',b1b2);
end

points = squirmer.getSwimmerBnd();
plot3(points(:,1),points(:,2),points(:,3))

shearFunc = @(t,x) shearFlow(x,0);
kernelPar = [1e-1,1];
treePar = {'nodeCapacity', 500, 'nearest', 1};

mobilityProb = @(t,x) mobilityProblem(squirmer, shearFunc,t,x,...
                                      kernelPar,treePar,20,0,1,1e-4,50,...
                                      [1e-5,100]);

ofunc = @(t,x) mobilityProb(t,x);                           
                                  
Totaltime = tic;
[t,odeout] = ode45(ofunc, [0 3], initialCon);
toc(Totaltime)

function vel = shearFlow(points, gradient)
vel = zeros(size(points));
vel(:,1) = gradient*points(:,3);
end

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