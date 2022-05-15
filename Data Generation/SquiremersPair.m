clear all; close all; clc;

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder))

pc = parcluster('local');
pc.JobStorageLocation = getenv('TMPDIR');
parpool(pc, str2num(getenv('SLURM_TASKS_PER_NODE')));

squirmer = Swimmers();
initialCon = [];

coarse = prolateSpheroid(0.5,0.5,0.2);
fine = prolateSpheroid(0.5,0.5,0.05);
fine = unique([fine;coarse],'rows','stable');
NN = nearestMatrix(coarse,fine,1,0);
vel = squirmerVel(coarse,2);

startlocal = [-1 0 0;
               1 0 0];

initialCon = [];
basis = Rotate(eye(3),0,7*pi/16,0,[0 0 0]);
b1b2 = basis(1:6)';
squirmer.addSwimmer(coarse,fine,NN,startlocal(1,:),b1b2,...
                    vel,[0;0;0],[0;0;0]);
initialCon = vertcat(initialCon,startlocal(1,:)',b1b2);

basis = Rotate(eye(3),0,-7*pi/16,0,[0 0 0]);
b1b2 = basis(1:6)';

squirmer.addSwimmer(coarse,fine,NN,startlocal(2,:),b1b2,...
                    vel,[0;0;0],[0;0;0]);
initialCon = vertcat(initialCon,startlocal(2,:)',b1b2);


[points,fine] = squirmer.getSwimmerBnd();
plot3(points(:,1),points(:,2),points(:,3))

shearFunc = @(t,x) shearFlow(x,0);
kernelPar = [1e-2,1];
treePar = {'nodeCapacity', 500, 'nearest', 1};

mobilityProb = @(t,x) mobilityProblemGravity(squirmer, shearFunc,t,x,...
                                      kernelPar,treePar,20,0,2,1e-4,100,...
                                      [1e-5,50],[0;0;-1;0.2]);

ofunc = @(t,x) mobilityProb(t,x);                           

filename = 'SquiremerNearPair.csv';
outfunc = @(t,y,flag) save(t,y,flag,filename);

opts = odeset('OutputFcn',outfunc);

Totaltime = tic;
[t,odeout] = ode45(ofunc, [0 2], initialCon,opts);
toc(Totaltime)

writematrix([t';odeout'],'SquiremerNearFullPair.csv')

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

function status = save(t,y,flag,filename)
    if ~isempty(flag)
        if flag == 'init'
            writematrix([],filename);
        end
    else
        data = readmatrix(filename);
        data = [data,[t;y]];
        writematrix(data,filename);
        disp(t)
    end
    status = 0;
end