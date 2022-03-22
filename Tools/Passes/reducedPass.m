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

potPoints = tree.points(tree.potentials,:);
targetPoints = tree.points(tree.targets,:);

nodes = (1:tree.nodeCount)';
nodes = nodes(isleaf(tree,nodes,'Target'))';

velpar = cell(numel(nodes),1);
vel = zeros(size(targetPoints));

parfor (i = 1:numel(nodes), arguments.parThreads)
    node = nodes(i);
    
    nodeSlice = tree.pointIndex(tree.targets) == node;
    leafpoints = targetPoints(nodeSlice,:);
    
    u = tree.interactions{node,1};
    
    if tree.arguments.nearest
        
        NN = sparse(0,0);
        usurf = zeros(0,3);
        upot = zeros(0,3);
        for ui = 1:size(u,2)
            unode = u(ui);

            uislice = tree.pointIndex == unode;

            upot = vertcat(upot,potentials(uislice,:));

            NNi = tree.NN(:,uislice);
            NNslice = any(NNi,2);
            NN = blkdiag(NN,NNi(NNslice',:));

            usurf = vertcat(usurf,tree.finePoints(NNslice',:));

        end

        leafpoints = reshape(leafpoints',[],1);
        NN = kron(NN,speye(3));
        usurf = reshape(usurf.',[],1);
        upot = reshape(upot',[],1);

        velpar{i} = DownleafNearest(leafpoints,zeros(0,1),zeros(0,1),...
                        usurf,upot,zeros(0,1),zeros(0,1),NN,blockSize,...
                        kernelPar,GPU);
    else

        usurf = zeros(tree.nodeCapacity*size(u,2),3);
        upot = zeros(tree.nodeCapacity*size(u,2),3);
        uindex = zeros(tree.nodeCapacity*size(u,2),1);
        for ui = 1:size(u,2)
            unode = u(ui);

            uislice = tree.pointIndex(tree.potentials) == unode;

            uipoints = potPoints(uislice,:);
            uipointpot = potentials(uislice,:);

            points = size(uipoints,1);

            usurf(((ui-1)*tree.nodeCapacity)+1:...
                  ((ui-1)*tree.nodeCapacity)+points,:) = ...
                  uipoints;
            upot(((ui-1)*tree.nodeCapacity)+1:...
                 ((ui-1)*tree.nodeCapacity)+points,:) = ...
                 uipointpot;
            uindex(((ui-1)*tree.nodeCapacity)+1:...
                 ((ui-1)*tree.nodeCapacity)+points,:) = ...
                 1;

        end
        usurf(~uindex,:) = [];
        upot(~uindex,:) = [];
        
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