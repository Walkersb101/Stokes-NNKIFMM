function [pot] = ReducedSolve(tree,vel,arguments,GPU)
%REDUCEDSOLVE Summary of this function goes here
%   Detailed explanation goes here

kernelPar = arguments.kernelPar;
blockSize = arguments.blockSize;

nodes = (1:tree.nodeCount)';
nodes = nodes(isleaf(tree,nodes,'Target'))';

blockStorage = cell(numel(nodes),1);
pot = zeros(sum(tree.targets),3);

potPoints = tree.points(tree.potentials,:);
targetPoints = tree.points(tree.targets,:);
finePoints = tree.finePoints;
index = tree.pointIndex(tree.targets);
NNfull = tree.NN;

NEAREST = tree.arguments.nearest;

if GPU
    
    velStorage = cell(numel(nodes),1);
    
    potPoints = gpuArray(potPoints);
    targetPoints = gpuArray(targetPoints);
    finePoints = gpuArray(finePoints);
    index = gpuArray(index);
    NNfull = gpuArray(NNfull);
    vel = gpuArray(vel);
    
    for i = 1:numel(nodes)
        node = nodes(i);

        nodeSlice = index == node;
        leafpoints = targetPoints(nodeSlice,:);

        if NEAREST

            uvel = vel(nodeSlice,:);

            NNi = NNfull(:,nodeSlice);
            NNslice = any(NNi,2);
            NN = NNi(NNslice',:);

            usurf = finePoints(NNslice',:);

            leafpoints = reshape(leafpoints',[],1);
            NN = kron(NN,speye(3));
            usurf = reshape(usurf.',[],1);
            uvel = reshape(uvel',[],1);

            blockStorage{i} = kernel(usurf,leafpoints,kernelPar)*NN;
            velStorage{i} = uvel;
        else

            usurf = potPoints(nodeSlice,:);
            uvel = vel(nodeSlice,:);

            leafpoints = reshape(leafpoints',[],1);
            usurf = reshape(usurf',[],1);
            uvel = reshape(uvel',[],1);
            
            blockStorage{i} = sparse(kernel(usurf,leafpoints,kernelPar));
            velStorage{i} = uvel;

        end
    end
    
    diagPot = gpuSpBlkDiag(blockStorage)\vertcat(velStorage{:});
    
    diagPot = reshape(diagPot,3,[])';
    
    count = 1;
    for i = 1:numel(nodes)
        node = nodes(i);

        nodeslice = tree.pointIndex(tree.targets) == node;

        pot(nodeslice,:) = gather(diagPot(count:count+sum(nodeslice)-1,:));
        count = count + sum(nodeslice);
    end
    
    
else
    
    parfor (i = 1:numel(nodes), arguments.parThreads)
        node = nodes(i);

        nodeSlice = index == node;
        leafpoints = targetPoints(nodeSlice,:);

        if NEAREST

            uvel = vel(nodeSlice,:);

            NNi = NNfull(:,nodeSlice);
            NNslice = any(NNi,2);
            NN = NNi(NNslice',:);

            usurf = finePoints(NNslice',:);

            leafpoints = reshape(leafpoints',[],1);
            NN = kron(NN,speye(3));
            usurf = reshape(usurf.',[],1);
            uvel = reshape(uvel',[],1);

            blockStorage{i} = reshape((kernel(usurf,leafpoints,...
                                              kernelPar)*NN)...
                                              \uvel,3,[])';
        else

            usurf = potPoints(nodeSlice,:);
            uvel = vel(nodeSlice,:);

            leafpoints = reshape(leafpoints',[],1);
            usurf = reshape(usurf',[],1);
            uvel = reshape(uvel',[],1);
            
            
            blockStorage{i} = reshape(kernel(usurf,leafpoints,kernelPar)...
                                      \uvel,3,[])';
        end
    end
    
    for i = 1:numel(nodes)
        node = nodes(i);

        nodeslice = tree.pointIndex(tree.targets) == node;

        pot(nodeslice,:) = blockStorage{i};
    end
end
end

% compute the index vector from each matrix.
function [p, q] = getIndexVectors(cellX)
numEl = numel(cellX);
p = zeros(1, numEl+1);
q = zeros(1, numEl+1);
for i = 1:numEl
    x = cellX{i};
    p(i+1) = size(x,1);
    q(i+1) = size(x,2);
end
p = cumsum(p);
q = cumsum(q);
end

function A = gpuSpBlkDiag(cellX)
    [p,q] = getIndexVectors(cellX);
    i=[];j=[];s=[];
    for n=1:numel(cellX)
        [in,jn,sn] = find(cellX{n});
        i = vertcat(i,in+p(n));
        j = vertcat(j,jn+q(n));
        s = vertcat(s,sn);
    end
    A = sparse(i,j,s);
end
