function [vel] = reducedPass(tree,potentials,arguments,GPU)
%reducedPass Compute downward pass of the KIFMM method on the reduced 
%   method.
%   The reduced preconditioning uses the same tree as regualar KIFMM
%   however only the interations between target and potential points within
%   the same node (and its children) are computed and all far field
%   interactions are ignored.
%   
% Inputs:
%   tree       : A tree struture containing potential points
%   potentials : A (N,3) array of potentials
%   uppot      : Uppwards equivent potentials calculated in the upward pass
%   arguments  : Struture passed from FMM class containing all variables, 
%                see FMM code for fields
% Output:
%   vel : A (M,3) array of velocities at the target points

kernelPar = arguments.kernelPar;
blockSize = arguments.blockSize;

nodes = (1:tree.nodeCount)';
nodes = nodes(isleaf(tree,nodes,'Target'))';

velpar = cell(numel(nodes),1);
vel = zeros(sum(tree.targets),3);

potPoints = tree.points(tree.potentials,:);
targetPoints = tree.points(tree.targets,:);
finePoints = tree.finePoints;
index = tree.pointIndex(tree.targets);
NNfull = tree.NN;

NEAREST = tree.arguments.nearest;

parfor (i = 1:numel(nodes), arguments.parThreads)
    node = nodes(i);
    
    nodeSlice = index == node;
    leafpoints = targetPoints(nodeSlice,:);
        
    if NEAREST
        
        upot = potentials(nodeSlice,:);

        NNi = NNfull(:,nodeSlice);
        NNslice = any(NNi,2);
        NN = NNi(NNslice',:);

        usurf = finePoints(NNslice',:);

        leafpoints = reshape(leafpoints',[],1);
        NN = kron(NN,speye(3));
        usurf = reshape(usurf.',[],1);
        upot = reshape(upot',[],1);

        velpar{i} = DownleafNearest(leafpoints,zeros(0,1),zeros(0,1),...
                        usurf,upot,zeros(0,1),zeros(0,1),NN,blockSize,...
                        kernelPar,GPU);
    else

        usurf = potPoints(nodeSlice,:);
        upot = potentials(nodeSlice,:);
        
        leafpoints = reshape(leafpoints',[],1);
        usurf = reshape(usurf',[],1);
        upot = reshape(upot',[],1);

        velpar{i} = Downleaf(leafpoints,zeros(0,1),zeros(0,1),usurf,...
                       upot,zeros(0,1),zeros(0,1),blockSize,kernelPar,GPU);

    end
end

for i = 1:numel(nodes)
    node = nodes(i);
    
    nodeslice = tree.pointIndex(tree.targets) == node;
    
    vel(nodeslice,:) = velpar{i};
end
end