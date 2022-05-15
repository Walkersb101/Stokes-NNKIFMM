clear all;

spmd
   gpuDevice(1); 
end

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder))

epsilon = 0.01;
mu = 1;

discretisation = 10:10:200;
samples = length(discretisation);

storage = zeros(8,samples);

for i = 1:samples
    disp(i+"\"+samples)
    
    N=discretisation(i);
    
    points = sphere(1,[0 0 0],N);    
    
    degreesoffreedom = numel(points);
    
    vel = ones(size(points,1),3);
    
    fmmtime = zeros(7,1);
    index = 0;
    for j = 0:2:12
        tic
        tree = OcTree(points,points,'nodeCapacity',500,'maxDepth',21);
        A = KIFMM(tree,[0.01,1],'GPU',1,'parThreads',j);

        A.computeVel(vel);
        fmmtime(index+1) = toc;
        index = index + 1;
        disp(fmmtime(index))
    end
    
    storage(:,i) = [numel(points); fmmtime];
    csvwrite('MathsGPUSphereBenchmark.csv',storage);
end
