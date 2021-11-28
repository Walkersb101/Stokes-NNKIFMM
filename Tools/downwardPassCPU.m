function [vel] = downwardPassCPU(tree,potentials,...
                                    uppot,arguments)
%downwardPassCPU Compute downward pass of the KIFMM method on the CPU 
%   The downward pass of the KIFMM method is comprised of two steps, first
%   is the computation of the downward equivent potentials from the
%   far field and parent nodes. This is achived through the computation of
%   the upwards potentials of far field nodes onto the downward surface and
%   then the calcuation of the downwards equivent potentials on the
%   downward surface. 
%   The for each leaf node the velocities of the target points are
%   calculated though the effect of the downwards equivent potentials and
%   point to point interactions of points in the near field
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

coronaRes = arguments.coronaRes;
coronaShells = arguments.coronaShells;
coronaPoints = coronaRes^3 - (coronaRes-2*coronaShells)^3;

levels = max(tree.nodeLevel);

downpot = zeros(tree.nodeCount,coronaPoints*3);

potPoints = tree.points(tree.potentials,:);
targetPoints = tree.points(tree.targets,:);

% compute from top down the downward equivalent potentials
for level = 1:1:levels
    nodes = find(tree.nodeLevel==level);

    downpottemp = zeros(size(nodes,2),coronaPoints*3);

   % compute downward potential of node of the level in parallel
    parfor (i = 1:size(nodes,2), arguments.parThreads)
        node = nodes(i);
        
        upsurf = genupsurf(tree,node,coronaRes,coronaShells);
        upsurf = reshape(upsurf.',[],1);

        downsurf = gendownsurf(tree,node,coronaRes,coronaShells);
        downsurf = reshape(downsurf.',[],1);
        
        v = tree.interactions{node,2};
        x = tree.interactions{node,4};
        
        RHS = zeros(coronaPoints*3,1);
        
        % compute effect upward potential of nodes in V
        if ~isempty(v)
            
            vsurf = zeros(coronaPoints*size(v,2),3);
            for vi = 1:size(v,2)
                vnode = v(vi);
                vsurf(((vi-1)*coronaPoints)+1:vi*coronaPoints,:) = ...
                    genupsurf(tree,vnode,coronaRes,coronaShells);
            end
            vsurf = reshape(vsurf',[],1);
            
            vpot = reshape(uppot(v,:)',[],1);
           
            RHS = RHS + blockcomputation(vsurf,upsurf,vpot,blockSize,...
                                         kernelPar);
            
        end
        
        % compute effect of potential points in X on downward equivent 
        % potential 
        if ~isempty(x)
            
            xsurf = zeros(tree.nodeCapacity*size(x,2),3);
            xpot = zeros(tree.nodeCapacity*size(x,2),3);
            for xi = 1:size(x,2)
                xnode = x(xi);

                xislice = tree.pointIndex(tree.potentials) == xnode;
                
                xipoints = potPoints(xislice,:);
                xipointpot = potentials(xislice,:);
                
                points = size(xipoints,1);
                
                xsurf(((xi-1)*tree.nodeCapacity)+1:...
                      ((xi-1)*tree.nodeCapacity)+points,:) = ...
                      xipoints;
                xpot(((xi-1)*tree.nodeCapacity)+1:...
                     ((xi-1)*tree.nodeCapacity)+points,:) = ...
                     xipointpot;
                
            end
            xsurf = xsurf(any(xsurf,2),:);
            xpot = xpot(any(xpot,2),:);
            
            xsurf = reshape(xsurf',[],1);
            xpot = reshape(xpot',[],1);
            
            RHS = RHS + blockcomputation(xsurf,upsurf,xpot,blockSize,...
                                         kernelPar);
            
        end
        
        % compute effect of parent on downward equivent potential
        parentnode = tree.nodeParents(node);
        
        parentsurf = gendownsurf(tree,parentnode,coronaRes,coronaShells);
        parentsurf = reshape(parentsurf',[],1);
        
        parentpot = downpot(parentnode,:)';
        
        RHS = RHS + blockcomputation(parentsurf,upsurf,parentpot,blockSize,kernelPar);
        
        % Compute potentials on downsurf from velocities on upward surfece 
        pot = kernel(downsurf,upsurf,kernelPar) \ RHS;

        downpottemp(i,:) = pot';
    end
    
    downpot(nodes,:) = downpottemp;
    
end

leaves = find(isleaf(tree,1:tree.nodeCount));

vel = zeros(size(potentials));
velpar = cell(size(leaves,1),1);

% compute potentials on points forall leaf nodes
parfor (i = 1:size(leaves,1), arguments.parThreads)
    
    node = leaves(i);
    
    downsurf = gendownsurf(tree,node,coronaRes,coronaShells);
    downsurf = reshape(downsurf',[],1);
    
    nodepot = downpot(node,:)';
    
    leafslice = tree.pointIndex(tree.targets) == node;
    leafpoints = targetPoints(leafslice,:);
    leafpoints = reshape(leafpoints',[],1);
    
    veltemp = blockcomputation(downsurf,leafpoints,nodepot,blockSize,kernelPar);
    
    u = tree.interactions{node,1};
    w = tree.interactions{node,3};
    
    % compute point to point iteraction for connected nodes
    if ~isempty(u)

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

        usurf = reshape(usurf',[],1);
        upot = reshape(upot',[],1);

        veltemp = veltemp + blockcomputation(usurf,leafpoints,upot,blockSize,kernelPar);

    end
    
    % compute effect of smaller nearby nodes
    if ~isempty(w)
            
        wsurf = zeros(coronaPoints*size(w,2),3);
        for wi = 1:size(w,2)
            wnode = w(wi);
            wsurf(((wi-1)*coronaPoints)+1:wi*coronaPoints,:) = ...
                genupsurf(tree,wnode,coronaRes,coronaShells);
        end
        wsurf = reshape(wsurf',[],1);

        wpot = reshape(uppot(w,:)',[],1);

        veltemp = veltemp + blockcomputation(wsurf,leafpoints,wpot,blockSize,kernelPar);

    end
    
    velpar{i} = reshape(veltemp,3,[])';
    
end

% map temporary variable back to final output 
for i = 1:size(leaves,1)
    
    node = leaves(i);
    
    leafslice = tree.pointIndex==node;

    vel(leafslice,:) = velpar{i};
    
end
end

