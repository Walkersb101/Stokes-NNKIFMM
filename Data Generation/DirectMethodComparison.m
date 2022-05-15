clear all;

pc = parcluster('local');
pc.JobStorageLocation = getenv('TMPDIR');
parpool(pc, str2num(getenv('SLURM_TASKS_PER_NODE')));

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder))

[X,Y] = meshgrid(15:10:115,1e-2);

storage = zeros(9,numel(X));

for i = 1:numel(X)
    
    disp( i + "/" + numel(X))
    
    points = sphere(1,[0 0 0],X(i));

    vel = [1 0 0].*ones(size(points,1),1);
    vel = reshape(vel.',[],1);
    
    tic
    tree = OcTree(points,points,'nodeCapacity',500);
    fmm = KIFMM(tree,[Y(i),1],'GPU',0,'parThreads',20,'GMRES',1);
    [pot1,~,~,iter1] = gmres(@(x) fmm.computeVel(x), vel,[],[],500);
    timeinverse1 = toc;
    
    pot1 = reshape(pot1,3,[]).';
    F = [sum(pot1(:,1));sum(pot1(:,2));sum(pot1(:,3))];
    
    err1 = norm(F - [6*pi; 0; 0])/norm([6*pi; 0; 0]);
    
    points = reshape(points.',[],1);
    
    tic
    [pot2,~,~,iter2] = gmres(@(x) blockcomputation(points,points,x,1,[Y(i),1]),vel,[],[],500);
    timeinverse2 = toc;

    pot2 = reshape(pot2,3,[]).';
    F = [sum(pot2(:,1));sum(pot2(:,2));sum(pot2(:,3))];
    err2 = norm(F - [6*pi; 0; 0])/norm([6*pi; 0; 0]);

    diff = norm(pot1-pot2)/norm(pot1);

    storage(:,i) = [numel(points); Y(i);timeinverse1;iter1(2);err1;...
                                        timeinverse2;iter2(2);err2;diff];
                                   
    clearvars points vel;
    csvwrite('DirectMethodComparison.csv',storage)
end