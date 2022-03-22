clear all; close all; clc;

spmd
   gpuDevice(1); 
end

trigrid=TriangularLattice(3,[10,10]);

spheres = Swimmers();
for i = 1:size(trigrid,1)
    spheres.addSwimmer(prolateSpheroid(1,2,[trigrid(i,1), trigrid(i,2), 0],0.1),[trigrid(i,1), trigrid(i,2), 0],[1;0;0;0;1;0]);
end

vel = [1 0 0];

N = sum(spheres.getSwimmerSizes());

ux = vel(1)*ones(N,1);
uy = vel(2)*ones(N,1);
uz = vel(3)*ones(N,1);

point_vel = interleaveComponents(ux,uy,uz);

points = spheres.getSwimmers();

spheres.genTree('nodeCapacity',2000);
spheres.genKIFMM([1e-5, 1],'parThreads',4,'GPU',1);

afunc = @(x) MobilityMatrixMulti(x, spheres);

out = gmres(afunc, [point_vel;zeros(6*spheres.swimmerNo(),1)],[],1e-4,100);

spheres.genKIFMM([1e-2, 1],'parThreads',4,'GPU',1,'GMRES',1);

Ax0 = afunc(out);
scale = (([point_vel;zeros(6*spheres.swimmerNo(),1)].')*Ax0)/norm(Ax0)^2;

out = gmres(afunc, [point_vel;zeros(6*spheres.swimmerNo(),1)],[],1e-4,100,[],[],scale*out);

