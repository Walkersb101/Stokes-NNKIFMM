clear all; close all; clc;

spmd
   gpuDevice(1); 
end

[X,Y] = TriangularLattice(0.4,[5,5]);
X = X - 2.5;
Y = Y - 2.5;

spheres = Swimmers();
initialCon = [];
for i = 1:numel(X)
    finePoints = prolateSpheroid(0.15,1,0.08);
    coarsePoints = finePoints(1:5:end,:);
    NN = nearestMatrix(reshape(coarsePoints.',[],1),reshape(finePoints.',[],1),2,4);
    spheres.addSwimmer(coarsePoints,finePoints,NN,[X(i), Y(i), 0],[1;0;0;0;1;0]);
    initialCon = vertcat(initialCon,[X(i), Y(i), 0]',[1;0;0;0;1;0]);
end

vel = [0 0 0];
mom = [0 0 0];

points = spheres.getSwimmerBnd();
Q = size(points,1);
x = points(:,1);
y = points(:,2);
z = points(:,3);

ux = vel(1)*ones(Q,1);
uy = vel(2)*ones(Q,1);
uz = vel(3)*ones(Q,1);

point_vel = [ux; uy; uz] + [mom(2) * z - mom(3) * y; mom(3) * x - mom(1) * z; mom(1) * y - mom(2) * x];
point_vel = interleaveComponents(point_vel(1:Q),point_vel(Q+1:2*Q),point_vel(2*Q+1:3*Q));

ofunc = @(t,x) SolveSystem(spheres, point_vel, t, x);

[t,y] = ode45(ofunc, [0 3], initialCon);

writematrix(vertcat(t',y'),"data.csv")

function [out] = SolveSystem(swimmers, vel, t, conditions)
out = zeros(size(conditions));
gradient = 0.2;

for i = 0:(numel(conditions)/9)-1
     data = conditions(i*9+1:(i+1)*9);
     swimmers.updateSwimmmer(i+1,data(1:3)',data(4:9))
end

swimmers.genTree('nodeCapacity',2000,'nearest',1);
swimmers.genKIFMM([1e-2, 1],'parThreads',4,'GPU',1,'GMRES',1,'blockSize',0.5);

shear = zeros(size(vel));
points = reshape((swimmers.Tree.points).',[],1);
shear(1:3:end) = gradient*points(3:3:end);

afunc = @(x) MobilityMatrixMulti(x, swimmers);

x = gmres(afunc, [vel-shear;zeros(6*swimmers.swimmerNo(),1)],[],1e-4,2000);

x0 = x(end-6*swimmers.swimmerNo()+1:end-3*swimmers.swimmerNo());
omega = x(end-3*swimmers.swimmerNo()+1:end);

for i = 0:swimmers.swimmerNo()-1
    out((i*9)+1:(i*9)+3) = x0((i*3)+1:(i*3)+3);
    out((i*9)+4:(i*9)+6) = cross(omega((i*3)+1:(i*3)+3),...
        swimmers.b{i+1,1});
    out((i*9)+7:(i*9)+9) = cross(omega((i*3)+1:(i*3)+3),...
        swimmers.b{i+1,2});
end
end