clear all;

pc = parcluster('local');
pc.JobStorageLocation = getenv('TMPDIR');
parpool(pc, str2num(getenv('SLURM_TASKS_PER_NODE')));

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder))

exact = [6*pi;0;0];

points = logspace(log10(0.04),0,10);
EPS = [1e-5 3e-04 1e-2 1e-1 0.5];


storage = zeros(11,numel(points)*numel(EPS));

for i = 1:numel(points)
    disp(i+"/"+numel(points))
    coarse = prolateSpheroid(1,1,points(i));
    fine = prolateSpheroid(1,1,points(i)/4);
    fine2 = [fine;coarse];
    
    vel = zeros(size(coarse));vel(:,1) = 1; vel = reshape(vel.',[],1);
    vel2 = zeros(size(fine));vel2(:,1) = 1; vel2 = reshape(vel2.',[],1);
    
    for j =1:numel(EPS)
        tic
        Tree = OcTree(coarse,coarse,'nodeCapacity',500);
        FMM = KIFMM(Tree,[EPS(j) 1],'parThreads',20,'blockSize',1,'GMRES',1);
        [F,~,~,iter1] = gmres(@(x) FMM.computeVel(x),vel,[],[],500);
        R = [sum(F(1:3:end));sum(F(2:3:end));sum(F(3:3:end))];
        errNys1 = norm(exact-R)/norm(exact);
        time1 = toc;

        tic
        NN = nearestMatrix(coarse,fine,2,4);
        Tree = OcTree(coarse,coarse,'nodeCapacity',500,'finePoints',fine,'nearest',1,'NN',NN);
        FMM = KIFMM(Tree,[EPS(j) 1],'parThreads',20,'blockSize',1,'GMRES',1);
        [F,~,~,iter3] = gmres(@(x) FMM.computeVel(x),vel,[],[],500);
        NN = kron(NN,speye(3));
        [Fx, Fy, Fz] = extractComponents(NN * F);
        R = [sum(Fx);sum(Fy);sum(Fz)];
        errNEAR1 = norm(exact-R)/norm(exact);
        time3 = toc;

        tic
        NN = nearestMatrix(coarse,fine2,2,4);
        Tree = OcTree(coarse,coarse,'nodeCapacity',500,'finePoints',fine2,'nearest',1,'NN',NN);
        FMM = KIFMM(Tree,[EPS(j) 1],'parThreads',20,'blockSize',1,'GMRES',1);
        [F,~,~,iter4] = gmres(@(x) FMM.computeVel(x),vel,[],[],500);
        NN = kron(NN,speye(3));
        [Fx, Fy, Fz] = extractComponents(NN * F);
        R = [sum(Fx);sum(Fy);sum(Fz)];
        errNEAR2 = norm(exact-R)/norm(exact);
        time4 = toc;
        
        data = [points(i),EPS(j),errNys1,errNEAR1,errNEAR2,time1,time3,time4,iter1(2),iter3(2),iter4(2)];
        disp(data)
        storage(:,(i-1)*numel(EPS) + j) = data.';
        writematrix(storage,'NearestComparionFMM.csv')
    end
end