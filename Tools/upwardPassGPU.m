function [uppot] = upwardPassGPU(tree,potentials,arguments)
%UPWARDPASSCPU Compute upwards pass of the KIFMM method on the GPU 
%   The upwards pass of the KIFMM method is compued onto the upwards
%   equivent surfaces though the calculation of the potential points onto
%   the downwards surface, where then the equivilent potentials are 
%   calculated on the upwards surface. The algorithm then works up
%   the tree calucating the equivent potential based on the equivelent
%   potentials of its children. This method uses gpuArray to passivly
%   acclerate the method using gpu computation.
% Inputs:
%   tree       : A tree struture containing potential points
%   potentials : A (N,3) array of potentials
%   arguments  : Struture passed from FMM class containing all variables, 
%                see FMM code for fields
% Output:
%   uppot : A (nodeCount,coronaPoints) array containg the equivent upwards
%           potential

kernelPar = arguments.kernelPar;
blockSize = arguments.blockSize;

coronaRes = arguments.coronaRes;
coronaShells = arguments.coronaShells;
coronaPoints = coronaRes^3 - (coronaRes-2*coronaShells)^3;

levels = max(tree.nodeLevel);
uppot = zeros(tree.nodeCount,coronaPoints*3);

potPoints = tree.points(tree.potentials,:);

% compute from bottom of the tree upwards
for level = levels:-1:0
    nodes = find(tree.nodeLevel==level);

    uppottemp = zeros(size(nodes,2),coronaPoints*3);

    % compute upward potential of node of the level in parallel
    for i = 1:size(nodes,2)
%     parfor (i = 1:size(nodes,2), arguments.parThreads)
        node = nodes(i);

        upsurf = genupsurf(tree,node,coronaRes,coronaShells);
        upsurf = gpuArray(reshape(upsurf.',[],1));

        downsurf = gpuArray(gendownsurf(tree,node,coronaRes,coronaShells));
        downsurf = reshape(downsurf.',[],1);

        % Compute potential points onto equivent surface for leaf nodes or
        % from the children equivalent surfaces 
        if isleaf(tree,node)

            nodeslice = tree.pointIndex(tree.potentials) == node;
            
            nodepoints = gpuArray(potPoints(nodeslice,:));
            nodepotentials = gpuArray(potentials(nodeslice,:));

            
            
            nodepoints = reshape(nodepoints.',[],1);
            nodepotentials = reshape(nodepotentials.',[],1);
            
            RHS = blockcomputation(nodepoints,downsurf,nodepotentials,...
                                   blockSize,kernelPar);

        else

            children = tree.nodeChildren(node,:);
            childrensurf = zeros(coronaPoints*8,3);
            for j = 1:8
                childrensurf(((j-1)*coronaPoints)+1:j*coronaPoints,:) = ...
                    gpuArray(genupsurf(tree,children(j),coronaRes,...
                             coronaShells));
            end
            childrensurf = reshape(childrensurf',[],1);

            childrenpot = gpuArray(uppot(children,:));
            childrenpot = reshape(childrenpot',[],1);

            RHS = blockcomputation(childrensurf,downsurf,childrenpot,...
                                   blockSize,kernelPar);

        end

        pot = kernel(upsurf,downsurf,kernelPar) \ RHS;

        uppottemp(i,:) = gather(pot');
    end

    uppot(nodes,:) = uppottemp;
end

