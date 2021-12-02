clear all;
% warning('off','MATLAB:nearlySingularMatrix');

epsilon = 0.01;
mu = 1;


discretisation = 47:3:65;
samples = length(discretisation);

storage = zeros(5,samples);

for i = 1:samples
    N=discretisation(i);
    
    gpuDevice(1);
    
    points = [sphere(1,[-10 -10 -10],N);
          sphere(1,[10 -10 -10],N);
          sphere(1,[-10 10 -10],N);
          sphere(1,[10 10 -10],N);
          sphere(1,[-10 -10 10],N);
          sphere(1,[10 -10 10],N);
          sphere(1,[-10 10 10],N);
          sphere(1,[10 10 10],N)];    
    
    disp(numel(points));
      
    potentials = 100*rand(size(points))-50;

    tic
    N = size(points,1);
    fine = gpuArray(reshape(points.',[],1));
    pointpotreshape = gpuArray(reshape(potentials.',[],1));

    M = size(fine,1)/3;

    blockNodes=floor(0.1*2^27/(9*M));
    blocks = [1:3*blockNodes:3*N (3*N)+1];

    directvel = zeros(3*N,1);
%     for j = 1:length(blocks)-1
%         directvel(blocks(j):blocks(j+1)-1,:) = kernel(fine,fine(blocks(j):blocks(j+1)-1,:),[0.01,1])*pointpotreshape;
%     end
    directvel = gather(reshape(directvel,3,[])');
    directtime = toc

    fine = [];
    pointpotreshape = [];

    tic
    tree = OcTree(points,points,'nodeCapacity',800,'maxDepth',21);
    

    A = KIFMM(tree,[0.01,1],'GPU',1,'parThreads',0,'coronaRes',8,'coronaShells',1);

    fmmvel = A.computeVel(potentials);
    fmmtime = toc

    error = norm(fmmvel - directvel,'fro')/norm(directvel,'fro');
    
    storage(:,i) = [numel(points); directtime; fmmtime; error; max(tree.nodeLevel)];
end

ax1 = subplot(2,1,1);
x = linspace(min(storage(1,1:end)),max(storage(1,:)),50);
plot(ax1,storage(1,1:end),storage(2,1:end),'b',storage(1,1:end),storage(3,1:end),'r',x,(1.5e-9)*x.^2,'k--',x,(1.4e-5)*x.*log(x),'k-.');

ax1.XScale = 'log';
ax1.YScale = 'log';

xlabel('Degrees of freedom') 
ylabel('Computational time [s]') 

xlim([min(storage(1,1:end))*0.9 max(storage(1,1:end))*1.1])

legend('Direct solver','KIFMM', 'N^2', 'NlogN','Location','northwest')

ax2 = subplot(2,1,2);
plot(storage(1,1:end),storage(4,1:end),'k')

ax2.XScale = 'log';
ax2.YScale = 'log';

xlabel('Degrees of freedom') 
ylabel('Relative Error between direct and KIFMM solution') 

xlim([min(storage(1,3:end))*0.9 max(storage(1,3:end))*1.1])
