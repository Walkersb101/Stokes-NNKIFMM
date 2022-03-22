% use to generate parallel Pool with 4 threads
% parpool(4);

% uncomment if using GPU
% spmd
%    gpuDevice(1); 
% end

clear all;


% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder))

% Generate Points on sphere in using cubic discretisation
points = sphere(1,[0 0 0],40);
fprintf("number of scalar DOF=%d\n",numel(points));

% generate velocity at each point
vel = [1 0 0].*ones(size(points,1),1);
vel = reshape(vel.',[],1); %reshape into correct format

% Generate adaptive Octree structure
% nodeCapacity = number of points in node
tree = OcTree(points,points,'nodeCapacity',2000);

% Generate FMM handles 
% First input is tree structure, then [epsilon, mu]
% Set 'GPU' = 1 for gpu Compute, ParThreads = number of parallel threads, I
% recomend 4-6 for GPU depending on the tree size
fmmInitial = KIFMM(tree,[1e-5,1],'GPU',0,'parThreads',20,'GMRES',1);
fmmActual = KIFMM(tree,[1e-2,1],'GPU',0,'parThreads',20,'GMRES',1);

tree = []; points =[]; % clear tree and points as now stored in KIFMM

% Direct to actual Epsilon Solution
[pot1,~,~,iter1] = gmres(@fmmActual.computeVel,vel,[],1e-6,50);

% Compute Inital guess
[pot21,~,~,iter21] = gmres(@fmmInitial.computeVel,vel,[],1e-6,50);

% Hegedus Trick
Ax0 = fmmActual.computeVel(pot21); % Compute vector matrix product
scale = ((vel.')*Ax0)/(norm(Ax0)^2);

% use Inital Guess
[pot2,~,~,iter22] = gmres(@fmmActual.computeVel,vel,[],1e-6,50,[],[],scale*pot21);


% Calculate error
pot1 = reshape(pot1,3,[]).';
F = [sum(pot1(:,1));sum(pot1(:,2));sum(pot1(:,3))];
err1 = norm(F - [6*pi; 0; 0])/norm([6*pi; 0; 0]);

pot2 = reshape(pot2,3,[]).';
F = [sum(pot2(:,1));sum(pot2(:,2));sum(pot2(:,3))];
err2 = norm(F - [6*pi; 0; 0])/norm([6*pi; 0; 0]);

fprintf("Direct to epsilon converged in %d iteration with error %f\n",...
       [iter1(2),err1]);
fprintf("Inital Guess converged in %d iterations followed by %d " + ...
"iterations to solve for actual epsilon,\n with a total error of %f\n",...
[iter21(2),iter22(2),err2]);