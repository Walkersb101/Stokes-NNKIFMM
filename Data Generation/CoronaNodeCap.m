clear all;

pc = parcluster('local');
pc.JobStorageLocation = getenv('TMPDIR');
parpool(pc, str2num(getenv('SLURM_TASKS_PER_NODE')));

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder))

epsilon = 0.01;
mu = 1;

discretisation = 10:15:100;
CoronaRes = [4 6 8];
nodeCapacity = [200 600 1000 1500 2000 3000];

[X, Y, Z] = meshgrid(discretisation,CoronaRes, nodeCapacity);
samples = numel(X);

storage = zeros(6,samples);

for i = 1:samples
    disp(i+"\"+samples)
    
    points = sphere(1,[0 0 0],X(i)); 
    
    vel = [1 0 0].*ones(size(points,1),1);
    vel = reshape(vel.',[],1);
    
    tic
    tree = OcTree(points,points,'nodeCapacity',Z(i),'maxDepth',21);
    A = KIFMM(tree,[0.01,1],'GPU',0,'parThreads',20,'coronaRes',Y(i),'coronaShells',1);

    [fmmpot,~,~,iter] = A.computePot(vel,[],[],500);
    fmmtime = toc;

    fmmpot = reshape(fmmpot,3,[]).';
    F = [sum(fmmpot(:,1));sum(fmmpot(:,2));sum(fmmpot(:,3))];

    fmmError = norm(F - [6*pi; 0; 0])/norm([6*pi; 0; 0]);
    fmmIterations = iter(2);
    
    storage(:,i) = [numel(points); Y(i); Z(i); fmmtime; fmmIterations; fmmError];
    csvwrite('CoronaNodeCap.csv',storage);
end