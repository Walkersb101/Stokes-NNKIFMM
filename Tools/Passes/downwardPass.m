function [vel] = downwardPass(tree,potentials,...
                                    uppot,arguments)
%downwardPass Compute downward pass of the KIFMM method on the CPU 
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
GPU = arguments.GPU;

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
            
        else
            vsurf = zeros(0,1);
            vpot = zeros(0,1);
        end
        
        % compute effect of parent on downward equivent potential
        parentnode = tree.nodeParents(node);
        
        parentsurf = gendownsurf(tree,parentnode,coronaRes,coronaShells);
        parentsurf = reshape(parentsurf',[],1);
        
        parentpot = downpot(parentnode,:)';
       
        if tree.arguments.nearest
            if ~isempty(x) 
                NN = sparse(0,0);
                xsurf = zeros(0,3);
                xpot = zeros(0,3);
                for xi = 1:size(x,2)
                    xnode = x(xi);

                    xislice = tree.pointIndex == xnode;
                
                    xpot = vertcat(xpot,potentials(xislice,:));
                    
                    NNi = tree.NN(:,xislice);
                    NNslice = any(NNi,2);
                    NN = blkdiag(sparse(NN),sparse(NNi(NNslice,:)));
                    
                    xsurf = vertcat(xsurf,tree.finePoints(NNslice',:));
                    
                end
                
                NN = kron(sparse(NN),speye(3));
                xsurf = reshape(xsurf.',[],1);
                xpot = reshape(xpot',[],1);

            else
                xsurf = zeros(0,1);
                xpot = zeros(0,1);
                NN = zeros(0);
            end

            downpottemp(i,:) = DownNonleafNearest(vsurf,vpot,xsurf,xpot,...
                parentsurf,parentpot,upsurf,downsurf,NN,blockSize,...
                kernelPar,GPU);
        else

            if ~isempty(x)
                
                xsurf = zeros(tree.nodeCapacity*size(x,2),3);
                xpot = zeros(tree.nodeCapacity*size(x,2),3);
                xslice = zeros(tree.nodeCapacity*size(x,2),1);
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
                    xslice(((xi-1)*tree.nodeCapacity)+1:...
                         ((xi-1)*tree.nodeCapacity)+points,:) = ...
                         1;

                end
                xslice = logical(xslice);
                xsurf = xsurf(xslice,:);
                xpot = xpot(xslice,:);

                xsurf = reshape(xsurf',[],1);
                xpot = reshape(xpot',[],1);

            else
                xsurf = zeros(0,1);
                xpot = zeros(0,1);
            end

            downpottemp(i,:) = DownNonleaf(vsurf,vpot,xsurf,xpot,...
                parentsurf,parentpot,upsurf,downsurf,blockSize,...
                kernelPar,GPU);
            
        end


    end
    
    downpot(nodes,:) = downpottemp;
    
end

leaves = find(isleaf(tree,1:tree.nodeCount,'Target'))';

vel = zeros(size(targetPoints));
velpar = cell(size(leaves,1),1);

% compute potentials on points for all leaf nodes
for i = 1:size(leaves,1)
%parfor (i = 1:size(leaves,1), arguments.parThreads)
    
    node = leaves(i);
    
    downsurf = gendownsurf(tree,node,coronaRes,coronaShells);
    downsurf = reshape(downsurf',[],1);
    
    nodepot = downpot(node,:)';
    
    leafslice = tree.pointIndex(tree.targets) == node;
    leafpoints = targetPoints(leafslice,:);
    leafpoints = reshape(leafpoints',[],1);
    
    u = tree.interactions{node,1};
    w = tree.interactions{node,3};
    
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
        
    else
        
        wsurf = zeros(0,1);
        wpot = zeros(0,1);
        
    end
    
    if tree.arguments.nearest
        if ~isempty(u) 
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

            NN = kron(NN,speye(3));
            usurf = reshape(usurf.',[],1);
            upot = reshape(upot',[],1);

        else
        
            usurf = zeros(0,1);
            upot = zeros(0,1);
            NN = zeros(0);
        
        end

        velpar{i} = DownleafNearest(leafpoints,downsurf,nodepot,...
                        usurf,upot,wsurf,wpot,NN,blockSize,kernelPar,GPU);
    else

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

        else
        
            usurf = zeros(0,1);
            upot = zeros(0,1);
        
        end

        velpar{i} = Downleaf(leafpoints,downsurf,nodepot,usurf,...
                        upot,wsurf,wpot,blockSize,kernelPar,GPU);

    end

end

% map temporary variable back to final output 
for i = 1:size(leaves,1)
    
    node = leaves(i);
    
    leafslice = tree.pointIndex(tree.targets) == node;

    vel(leafslice,:) = velpar{i};
    
end
end

