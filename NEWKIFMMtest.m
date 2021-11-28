clear all;
% warning('off','MATLAB:nearlySingularMatrix');

epsilon = 0.01;
mu = 1;

points = [sphere(1,[-10 -10 -10],30);
          sphere(1,[10 -10 -10],30);
          sphere(1,[-10 10 -10],30);
          sphere(1,[10 10 -10],30);
          sphere(1,[-10 -10 10],30);
          sphere(1,[10 -10 10],30);
          sphere(1,[-10 10 10],30);
          sphere(1,[10 10 10],30)];    
      
     
potentials = 100*rand(size(points))-50;

tic
N = size(points,1);
fine = gpuArray(reshape(points.',[],1));
pointpotreshape = gpuArray(reshape(potentials.',[],1));

Max_array_size = (0.1 * 8000000000)/64;
Row_block = floor(Max_array_size/(N*3));
blocks = [1:3*Row_block:3*N (3*N)+1];

directvel = zeros(3*N,1);
for i = 1:length(blocks)-1
    directvel(blocks(i):blocks(i+1)-1,:) = kernel(fine,fine(blocks(i):blocks(i+1)-1,:),[0.01,1])*pointpotreshape;
end
directvel = gather(reshape(directvel,3,[])');
toc

fine = [];
pointpotreshape = [];

disp("Nystrom done!")

tic
tree = OcTree(points,points,'nodeCapacity',800,'maxDepth',21);

A = KIFMM(tree,[0.01,1],'GPU',1,'parThreads',2);

fmmvel = A.computeVel(potentials);
toc

diff = fmmvel - directvel;
norm(diff,'fro')/norm(directvel,'fro')