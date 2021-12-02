clear all;

exact = blkdiag(6*pi*1*1*eye(3),8*pi*1*1*eye(3));

N = 30;
points = sphere(1,[0 0 0],N);

numel(points)

tree = OcTree(points,points,'nodeCapacity',800);

fmm = KIFMM(tree,[0.01,1],'GPU',1,'parThreads',2,'GMRES',1);

tic
R = Gen_Resistance(fmm);
toc

err = norm(exact-R,'fro')/norm(exact,'fro');


