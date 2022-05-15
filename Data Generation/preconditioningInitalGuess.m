clear all; close all; clc;

pc = parcluster('local');
pc.JobStorageLocation = getenv('TMPDIR');
parpool(pc, str2num(getenv('SLURM_TASKS_PER_NODE')));

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder))

[X,Y] = triangularLattice(0.4,[0,0],[5.2,5.1]);

spheres = Swimmers();
initialCon = [];
for i = 1:numel(X)
    points = prolateSpheroid(0.15,1,0.022);
    vel = zeros(size(points));
    spheres.addSwimmer(points,[X(i), Y(i), 0],[1;0;0;0;1;0],vel,...
                       [0;0;0],[0;0;0]);
    initialCon = vertcat(initialCon,[X(i), Y(i), 0]',[1;0;0;0;1;0]);
end
points = spheres.getSwimmerBnd();
N = size(points,1);
point_vel = zeros(3*N,1);

% data = SolveSystem(spheres, point_vel, 0, initialCon);
% disp("Non Pre done")
dataPre = SolveSystemPre(spheres, point_vel, 0, initialCon);
disp("Pre done")

% file = zeros(2,max(numel(data),numel(dataPre)));
% file(1,1:numel(data)) = data.';file(2,1:numel(dataPre)) = dataPre.';

writematrix(dataPre.',"PreconSmallInital.csv")

% [X,Y] = triangularLattice(0.4,[0,0],[16.6,16.6]);
% 
% spheres = Swimmers();
% initialCon = [];
% for i = 1:numel(X)
%     points = prolateSpheroid(0.15,1,0.18);
%     vel = zeros(size(points));
%     spheres.addSwimmer(points,[X(i), Y(i), 0],[1;0;0;0;1;0],vel,...
%                        [0;0;0],[0;0;0]);
%     initialCon = vertcat(initialCon,[X(i), Y(i), 0]',[1;0;0;0;1;0]);
% end
% points = spheres.getSwimmerBnd();
% N = size(points,1);
% point_vel = zeros(3*N,1);
% 
% data = SolveSystem(spheres, point_vel, 0, initialCon);
% disp("Non Pre done")
% dataPre = SolveSystemPre(spheres, point_vel, 0, initialCon);
% disp("Pre done")
% 
% file = zeros(2,max(numel(data),numel(dataPre)));
% file(1,1:numel(data)) = data.';file(2,1:numel(dataPre)) = dataPre.';
% 
% writematrix(file,"PreconLargeInital.csv")

quit

function [out] = SolveSystem(swimmers, vel, t, conditions)
tic
gradient = 2;

for i = 0:(numel(conditions)/9)-1
     data = conditions(i*9+1:(i+1)*9);
     swimmers.updateSwimmmer(i+1,data(1:3)',data(4:9))
end

swimmers.genTree('nodeCapacity',500,'nearest',0);
swimmers.genKIFMM([1e-5, 1],'parThreads',20,'GPU',0,'GMRES',1,'format',2,'blockSize',1);

shear = zeros(size(vel));
points = reshape((swimmers.Tree.points),[],1);
N=numel(points)/3;
shear(1:N) = gradient*points(2*N+1:end);

[B,BT] = genMatrixComponets(swimmers);

afunc = @(x) mobilityMatrixMulti(x,swimmers,B,BT);

input = [vel-shear;zeros(6*swimmers.swimmerNo(),1)];

[xp,~,~,~,out1] = gmres(afunc, input, [], 1e-6, 1000);
out1(end+1) = 0;

swimmers.genKIFMM([1e-2, 1],'parThreads',20,'GPU',0,'GMRES',1,'format',2,'blockSize',1);

Axp = afunc(xp);
scale = (input.'*Axp)/(norm(Axp)^2);

[~,~,~,~,out2] = gmres(afunc, input, [], 1e-6, 1000,[],[],scale*xp);


out = [out1;out2]/norm(input);
time = toc;
out = [time;out];
end

function [out] = SolveSystemPre(swimmers, vel, t, conditions)
tic
gradient = 2;

for i = 0:(numel(conditions)/9)-1
     data = conditions(i*9+1:(i+1)*9);
     swimmers.updateSwimmmer(i+1,data(1:3)',data(4:9))
end

swimmers.genTree('nodeCapacity',500,'nearest',0);
swimmers.genKIFMM([1e-5, 1],'parThreads',20,'GPU',0,'GMRES',1,'format',2,'blockSize',1);

shear = zeros(size(vel));
points = reshape((swimmers.Tree.points),[],1);
N=numel(points)/3;
shear(1:N) = gradient*points(2*N+1:end);

[B,BT,~] = genMatrixComponets(swimmers);

afunc = @(x) mobilityMatrixMulti(x,swimmers,B,BT);

input = [vel-shear;zeros(6*swimmers.swimmerNo(),1)];

[xp,~,~,~,out1] = gmres(afunc, input, [], 1e-6, 1000);
out1(end+1) = 0;

swimmers.genKIFMM([1e-2, 1],'parThreads',20,'GPU',0,'GMRES',1,'format',2,'blockSize',1);

xp(1:3*N) =  (1e-5/1e-2)*xp(1:3*N);

Axp = afunc(xp);
scale = (input.'*Axp)/(norm(Axp)^2);

[~,~,~,~,out2] = gmres(afunc, input, [], 1e-6, 1000,[],[],scale*xp);

out = [out1;out2]/norm(input);

time = toc;
out = [time;out];
end