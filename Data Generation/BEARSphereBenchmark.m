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

discretisation = 10:10:200;
samples = length(discretisation);

storage = zeros(12,samples);

for i = 1:15
    disp(i+"\"+samples)
    
    N=discretisation(i);
    
    points = sphere(1,[0 0 0],N);    
    
    degreesoffreedom = numel(points);
    
    vel = ones(size(points,1),3);
    
    fmmtime = zeros(11,1);
    index = 0;
    for j = 0:2:20
        tic
        tree = OcTree(points,points,'nodeCapacity',500,'maxDepth',21);
        A = KIFMM(tree,[0.01,1],'GPU',0,'parThreads',j,'coronaRes',6,'coronaShells',1);

        A.computeVel(vel);
        fmmtime(index+1) = toc;
        index = index + 1;
    end
    
    storage(:,i) = [numel(points); fmmtime];
    csvwrite('BEARSphereBenchmark.csv',storage);
end

