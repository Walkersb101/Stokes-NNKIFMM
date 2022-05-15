clear all; close all; clc;

pc = parcluster('local');
pc.JobStorageLocation = getenv('TMPDIR');
parpool(pc, str2num(getenv('SLURM_TASKS_PER_NODE')));

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder))

squirmer = Swimmers();
initialCon = [];

bndCoarse = [wallTri([0 0 0],[10 10],3,0.25);...
             wallTri([0 0 10],[10 10],3,0.25)];
bndFine = [wallTri([0 0 0],[10 10],3,0.1);...
             wallTri([0 0 10],[10 10],3,0.1)];
               
bndFine = unique([bndFine;bndCoarse],'rows','stable');
               
NN = nearestMatrix(bndCoarse,bndFine,4,5);
squirmer.updateBnd(bndCoarse,bndFine,NN,zeros(size(bndCoarse)));
disp([numel(bndCoarse), numel(bndFine)])
% bndFine = []; bndCoarse = [];


coarse = prolateSpheroid(0.5,0.5,0.33);
fine = prolateSpheroid(0.5,0.5,0.1);
fine = unique([fine;coarse],'rows','stable');
NN = nearestMatrix(coarse,fine,1,0);
vel = squirmerVel(coarse,1.5);

NoSquirmers = 50;

startlocal = [-1.65,-1.65,1] + [3.3,3.3,5].*rand(NoSquirmers,3);

initialCon = [];
for i = 1:NoSquirmers
    basis = Rotate(eye(3),2*pi*rand(1),2*pi*rand(1),2*pi*rand(1),[0 0 0]);
    b1b2 = basis(1:6)';
    
    squirmer.addSwimmer(coarse,fine,NN,startlocal(i,:),b1b2,...
                        vel,[0;0;0],[0;0;0]);
    initialCon = vertcat(initialCon,startlocal(i,:)',b1b2);
end

[points,fine,~] = squirmer.getSwimmerBnd();
disp([numel(points), numel(fine)])

shearFunc = @(t,x) shearFlow(x,0);
kernelPar = [1e-2,1];
treePar = {'nodeCapacity', 500, 'nearest', 1};

mobilityProb = @(t,x) mobilityProblemGravity(squirmer, shearFunc,t,x,...
                                      kernelPar,treePar,20,0,1,1e-4,75,...
                                      [2e-5,100],[0;0;-1;0.5]);

ofunc = @(t,x) mobilityProb(t,x);                           

filename = 'SquiremerGrav.csv';
outfunc = @(t,y,flag) save(t,y,flag,filename);

opts = odeset('OutputFcn',outfunc,'relTol',1e-3);

Totaltime = tic;
[t,odeout] = ode45(ofunc, [0 10], initialCon,opts);
toc(Totaltime)

writematrix([t';odeout'],'SquiremerGravFull.csv')

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