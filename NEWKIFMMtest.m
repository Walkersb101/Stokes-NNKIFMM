clear all;
% warning('off','MATLAB:nearlySingularMatrix');

epsilon = 0.01;
mu = 1;

points = [Sphere(1,[-10 -10 -10],10);
          Sphere(1,[10 -10 -10],10);
          Sphere(1,[-10 10 -10],10);
          Sphere(1,[10 10 -10],10);
          Sphere(1,[-10 -10 10],10);
          Sphere(1,[10 -10 10],10);
          Sphere(1,[-10 10 10],10);
          Sphere(1,[10 10 10],10)];    
      
     
potentials = 10*rand(size(points))-5;

disp(numel(points))


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
directvel = reshape(directvel,3,[])';
toc

disp("Nystrom done!")

tree = OcTree(points,'nodeCapacity',800,'maxDepth',21);

A = KIFMM(tree,[0.01,1]);

fmmvel = A.computeVel(potentials);

diff = fmmvel - directvel;
norm(diff,'fro')/norm(directvel,'fro')