function [uppot] = upwardPass(tree,potentials,arguments)
%UPWARDPASS Compute upwards pass of the KIFMM method on the CPU 
%   The upwards pass of the KIFMM method is compued onto the upwards
%   equivent surfaces though the calculation of the potential points onto
%   the downwards surface, where then the equivilent potentials are 
%   calculated on the upwards surface. The algorithm then works up
%   the tree calucating the equivent potential based on the equivelent
%   potentials of its children. 
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
GPU = arguments.GPU;

coronaRes = arguments.coronaRes;
coronaShells = arguments.coronaShells;
coronaPoints = size(gendownsurf(tree,1,coronaRes,coronaShells),1);

levels = max(tree.nodeLevel);
uppot = zeros(tree.nodeCount,coronaPoints*3);

potPoints = tree.points(tree.potentials,:);

% compute from bottom of the tree upwards
for level = levels:-1:0
    nodes = find(tree.nodeLevel==level & ...
                (isleaf(tree,1:tree.nodeCount,'Potential') | ...
                ~isleaf(tree,(1:tree.nodeCount)','All')'));
    
    uppottemp = zeros(size(nodes,2),coronaPoints*3);
           
    % compute upward potential of node of the level in parallel
    parfor (i = 1:size(nodes,2), arguments.parThreads)
    %for i = 1:size(nodes,2)
        node = nodes(i);

        upsurf = genupsurf(tree,node,coronaRes,coronaShells);
        upsurf = reshape(upsurf.',[],1);

        downsurf = gendownsurf(tree,node,coronaRes,coronaShells);
        downsurf = reshape(downsurf.',[],1);
        
        % Compute potential points onto equivent surface for leaf nodes or
        % from the children equivalent surfaces 
        if isleaf(tree,node,'Potential')

            nodeslice = tree.pointIndex(tree.potentials) == node;
            nodePotentials = potentials(nodeslice,:);
            nodePotentials = reshape(nodePotentials.',[],1);
            
            if tree.arguments.nearest
                
                NN = tree.NN(:,nodeslice);
                NNslice = any(NN,2);
                NN = NN(NNslice,:);
                NN = kron(sparse(NN),speye(3));

                nodePoints = tree.finePoints(NNslice',:);
                nodePoints = reshape(nodePoints.',[],1);
                
                uppottemp(i,:) = UpNearest(nodePoints,upsurf,downsurf,...
                          nodePotentials,NN,blockSize,kernelPar,GPU);
            else
                nodePoints = potPoints(nodeslice,:);
                nodePoints = reshape(nodePoints.',[],1);
                
                uppottemp(i,:) = Up(nodePoints,upsurf,downsurf,...
                                nodePotentials,blockSize,kernelPar,GPU);
            
            end

        else
            
            children = tree.nodeChildren(node,:);
            childrensurf = zeros(coronaPoints*8,3);
            for j = 1:8
                childrensurf(((j-1)*coronaPoints)+1:j*coronaPoints,:) = ...
                    genupsurf(tree,children(j),coronaRes,coronaShells);
            end
            childrensurf = reshape(childrensurf',[],1);

            childrenpot = uppot(children,:);
            childrenpot = reshape(childrenpot',[],1);
                               
            uppottemp(i,:) = Up(childrensurf,upsurf,downsurf,...
                                childrenpot,blockSize,kernelPar,GPU);

        end
    end
    
    uppot(nodes,:) = uppottemp;
    
end