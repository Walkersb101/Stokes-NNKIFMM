clear all; close all; clc;

spmd
   gpuDevice(1); 
end

[X,Y] = triangularLattice(0.4,[5,5],[0,0]);
X = X - 2.5;
Y = Y - 2.5;

spheres = Swimmers();
initialCon = [];
for i = 1:numel(X)
    finePoints = prolateSpheroid(0.15,1,0.08);
    coarsePoints = finePoints(1:4:end,:);
    NN = nearestMatrix(coarsePoints, finePoints, 1, 20);
    
    vel = zeros(size(finePoints));
    
    spheres.addSwimmer(finePoints,[X(i),Y(i),0],[1;0;0;0;1;0],...,
                       vel,[0;0;0],[0;0;0]);
    initialCon = vertcat(initialCon,[X(i), Y(i), 0]',[1;0;0;0;1;0]);
end

shearFunc = @(t,x) shearFlow(x,1);
kernelPar = [1e-2,1];
treePar = {'nodeCapacity', 500, 'nearest', 0};

mobilityProb = @(t,x) mobilityProblem(spheres, shearFunc,t,x,...
                                      kernelPar,treePar,4,1,1e-4,50,...
                                      [1e-6,100]);

Totaltime = tic;
mobilityProb(0,initialCon);
toc(Totaltime)

function vel = shearFlow(points, gradient)
vel = zeros(size(points));
vel(:,1) = gradient*points(:,3);
end

