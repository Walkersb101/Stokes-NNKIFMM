function [uppot] = upwardPassCPU(tree,potentials,arguments)
%UPWARDPASSCPU Compute upwards pass of the KIFMM method on the CPU 
%   The upwards pass of the KIFMM method is compued onto the upwards
%   equivent surfaces though the calculation of the potential nodes onto
%   the downwards check surface, where then the equivilent potentials are 
%   calculated on the upwards equivent points. The algorithm then works up
%   the tree calucating the equivent potential based on the equivelent
%   potentials of its children. 
% Inputs:
%   tree       : A tree struture containing potential points
%   potentials : 
%   arguments  : 

kernalPar = arguments.kernalPar;

coronaRes = arguments.coronaRes;
coronaShells = arguments.coronaShells;
coronapoints = coronaRes^3 - (coronaRes-2*coronaShells)^3;

levels = max(tree.nodeLevel);
uppot = zeros(tree.nodeCount,coronapoints*3);

for level = levels:-1:0
    nodes = find(tree.nodeLevel==level);

    uppottemp = zeros(size(nodes,2),coronapoints*3);

    parfor (i = 1:size(nodes,2), arguments.parThreads)
        node = nodes(i);

        upsurf = genupsurf(tree,node,coronaRes,coronaShells);
        upsurf = reshape(upsurf.',[],1);

        downsurf = gendownsurf(tree,node,coronaRes,coronaShells);
        downsurf = reshape(downsurf.',[],1);

        if isleaf(tree,node)

            nodeslice = tree.pointIndex == node;
            nodepoints = tree.points(nodeslice,:);
            nodepotentials = potentials(nodeslice,:);

            nodepoints = reshape(nodepoints.',[],1);
            nodepotentials = reshape(nodepotentials.',[],1);

            RHS = blockcomputation(nodepoints,downsurf,nodepotentials,kernalPar);

        else

            children = tree.nodeChildren(node,:);
            childrensurf = zeros(coronapoints*8,3);
            for j = 1:8
                childrensurf(((j-1)*coronapoints)+1:j*coronapoints,:) = ...
                    genupsurf(tree,children(j),coronaRes,coronaShells);
            end
            childrensurf = reshape(childrensurf',[],1);

            childrenpot = uppot(children,:);
            childrenpot = reshape(childrenpot',[],1);

            RHS = blockcomputation(childrensurf,downsurf,childrenpot,kernalParam);

        end

        pot = kernel(upsurf,downsurf,kernalParam) \ RHS;

        uppottemp(i,:) = pot';
    end

    uppot(nodes,:) = uppottemp;
end

