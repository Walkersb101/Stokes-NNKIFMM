function [vel] = downwardPassGPU(tree,potentials,...
                                    uppot,arguments)
%downwardPassCPU Compute downward pass of the KIFMM method on the GPU 
%   The downward pass of the KIFMM method is comprised of two steps, first
%   is the computation of the downward equivent potentials from the
%   far field and parent nodes. This is achived through the computation of
%   the upwards potentials of far field nodes onto the downward surface and
%   then the calcuation of the downwards equivent potentials on the
%   downward surface. 
%   The for each leaf node the velocities of the target points are
%   calculated though the effect of the downwards equivent potentials and
%   point to point interactions of points in the near field.
%   This method uses gpuArray to passivly accelrate the method using gpu 
%   computation.
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
        
        upsurf = gpuArray(genupsurf(tree,node,coronaRes,coronaShells));
        upsurf = reshape(upsurf.',[],1);

        downsurf = gpuArray(gendownsurf(tree,node,coronaRes,coronaShells));
        downsurf = reshape(downsurf.',[],1);
        
        v = tree.interactions{node,2};
        x = tree.interactions{node,4};
        
        RHS = zeros(coronaPoints*3,1,'gpuArray');
        
        % compute effect upward potential of nodes in V
        if ~isempty(v)
            
            vsurf = zeros(coronaPoints*size(v,2),3,'gpuArray');
            for vi = 1:size(v,2)
                vnode = v(vi);
                vsurf(((vi-1)*coronaPoints)+1:vi*coronaPoints,:) = ...
                    gpuArray(genupsurf(tree,vnode,coronaRes,coronaShells));
            end
            vsurf = reshape(vsurf',[],1);
            
            vpot = reshape(gpuArray(uppot(v,:)'),[],1);
           
            RHS = RHS + blockcomputation(vsurf,upsurf,vpot,blockSize,...
                                         kernelPar);
            
            vsurf = [];
            vpot = [];
        end
        
        % compute effect of potential points in X on downward equivent 
        % potential 
        if ~isempty(x)
            
            xsurf = zeros(tree.nodeCapacity*size(x,2),3,'gpuArray');
            xpot = zeros(tree.nodeCapacity*size(x,2),3,'gpuArray');
            for xi = 1:size(x,2)
                xnode = x(xi);

                xislice = tree.pointIndex(tree.potentials) == xnode;
                
                xipoints = gpuArray(potPoints(xislice,:));
                xipointpot = gpuArray(potentials(xislice,:));
                
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
            
                                     
            xpot = [];
            xsurf = [];
        end
        
        % compute effect of parent on downward equivent potential
        parentnode = tree.nodeParents(node);
        
        parentsurf = gpuArray(gendownsurf(tree,parentnode,coronaRes,...
                              coronaShells));
        parentsurf = reshape(parentsurf',[],1);
        
        parentpot = gpuArray(downpot(parentnode,:)');
        
        RHS = RHS + blockcomputation(parentsurf,upsurf,parentpot,...
                                     blockSize,kernelPar);
        
        % Compute potentials on downsurf from velocities on upward surfece 
        pot = kernel(downsurf,upsurf,kernelPar) \ RHS;

        downpottemp(i,:) = gather(pot');
    end
    
    downpot(nodes,:) = downpottemp;
    
end

leaves = find(isleaf(tree,1:tree.nodeCount));

vel = zeros(size(potentials));
velpar = cell(size(leaves,1),1);

% compute potentials on points forall leaf nodes
parfor (i = 1:size(leaves,1), arguments.parThreads)
    
    node = leaves(i);
    
    downsurf = gpuArray(gendownsurf(tree,node,coronaRes,coronaShells));
    downsurf = reshape(downsurf',[],1);
    
    nodepot = gpuArray(downpot(node,:)');
    
    leafslice = tree.pointIndex(tree.targets) == node;
    leafpoints = gpuArray(targetPoints(leafslice,:));
    leafpoints = gpuArray(reshape(leafpoints',[],1));
    
    veltemp = blockcomputation(downsurf,leafpoints,nodepot,blockSize,...
                               kernelPar);
    
    u = tree.interactions{node,1};
    w = tree.interactions{node,3};
    
    % compute point to point iteraction for connected nodes
    if ~isempty(u)

        usurf = zeros(tree.nodeCapacity*size(u,2),3,'gpuArray');
        upot = zeros(tree.nodeCapacity*size(u,2),3,'gpuArray');
        for ui = 1:size(u,2)
            unode = u(ui);

            uislice = tree.pointIndex(tree.potentials) == unode;

            uipoints = gpuArray(potPoints(uislice,:));
            uipointpot = gpuArray(potentials(uislice,:));

            points = size(uipoints,1);

            usurf(((ui-1)*tree.nodeCapacity)+1:...
                  ((ui-1)*tree.nodeCapacity)+points,:) = ...
                  uipoints;
            upot(((ui-1)*tree.nodeCapacity)+1:...
                 ((ui-1)*tree.nodeCapacity)+points,:) = ...
                 uipointpot;

        end
        
        uipoints = [];
        uipointpot = [];
        
        usurf(~any(usurf,2),:) = [];
        upot(~any(upot,2),:) = [];

        usurf = reshape(usurf',[],1);
        upot = reshape(upot',[],1);

        veltemp = veltemp + blockcomputation(usurf,leafpoints,upot,...
                                             blockSize,kernelPar);
                                         
        usurf = [];
        upot = [];

    end
    
    % compute effect of smaller nearby nodes
    if ~isempty(w)
            
        wsurf = zeros(coronaPoints*size(w,2),3,'gpuArray');
        for wi = 1:size(w,2)
            wnode = w(wi);
            wsurf(((wi-1)*coronaPoints)+1:wi*coronaPoints,:) = ...
                gpuArray(genupsurf(tree,wnode,coronaRes,coronaShells));
        end
        wsurf = reshape(wsurf',[],1);

        wpot = reshape(gpuArray(uppot(w,:)'),[],1);

        veltemp = veltemp + blockcomputation(wsurf,leafpoints,wpot,...
                                             blockSize,kernelPar);

        wsurf = [];
        wpot = [];
    end
    
    velpar{i} = gather(reshape(veltemp,3,[])');
    
end

% map temporary variable back to final output 
for i = 1:size(leaves,1)
    
    node = leaves(i);
    
    leafslice = tree.pointIndex==node;

    vel(leafslice,:) = velpar{i};
    
end
end

