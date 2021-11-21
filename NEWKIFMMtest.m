clear all;
% warning('off','MATLAB:nearlySingularMatrix');

epsilon = 0.01;
mu = 1;

points = [Sphere(1,[-10 -10 -10],10);
          Sphere(1,[10 -10 -10],10);
          Sphere(1,[-10 10 -10],10);
          Sphere(1,[10 10 -10],10);
          Sphere(1,[-10 -10 10],10);
          Sphere(1,[10 -10 10],10);
          Sphere(1,[-10 10 10],10);
          Sphere(1,[10 10 10],10)];    
      
     
potentials = 10*rand(size(points))-5;

disp(numel(points))


tic
N = size(points,1);
fine = gpuArray(reshape(points.',[],1));
pointpotreshape = gpuArray(reshape(potentials.',[],1));

Max_array_size = (0.1 * 8000000000)/64;
Row_block = floor(Max_array_size/(N*3));
blocks = [1:3*Row_block:3*N (3*N)+1];

nystromvel = zeros(3*N,1);
for i = 1:length(blocks)-1
    nystromvel(blocks(i):blocks(i+1)-1,:) = s(fine,fine(blocks(i):blocks(i+1)-1,:),epsilon,mu)*pointpotreshape;
end
nystromvel = reshape(nystromvel,3,[])';
toc

disp("Nystrom done!")


tic

coronares = 6;
coronarings = 1;
coronapoints = coronares^3 - (coronares-2*coronarings)^3;

tree = OcTree(points,'nodeCapacity',400,'maxDepth',21);

points = [];
toc

levels = max(tree.nodeLevel);
uppot = zeros(tree.nodeCount,coronapoints*3);
downpot = zeros(tree.nodeCount,coronapoints*3);


tic
for level = levels:-1:0
    nodes = find(tree.nodeLevel==level);
    
    uppottemp = zeros(size(nodes,2),coronapoints*3);
    
    for i = 1:size(nodes,2)
        node = nodes(i);
            
        upsurf = genupsurf(tree,node,coronares,coronarings);
        upsurf = reshape(upsurf.',[],1);

        downsurf = gendownsurf(tree,node,coronares,coronarings);
        downsurf = reshape(downsurf.',[],1);
        
        if isleaf(tree,node)
            
            nodeslice = tree.pointIndex == node;
            nodepoints = tree.points(nodeslice,:);
            nodepotentials = potentials(nodeslice,:);
            
            nodepoints = reshape(nodepoints.',[],1);
            nodepotentials = reshape(nodepotentials.',[],1);
            
%             RHS = s(nodepoints,downsurf,epsilon,mu)*nodepotentials;
            RHS = blockcomputation(nodepoints,downsurf,epsilon,mu,nodepotentials);
        
        else
            
            children = tree.nodeChildren(node,:);
            childrensurf = zeros(coronapoints*8,3);
            for j = 1:8
                childrensurf(((j-1)*coronapoints)+1:j*coronapoints,:) = ...
                    genupsurf(tree,children(j),coronares,coronarings);
            end
            childrensurf = reshape(childrensurf',[],1);
            
            childrenpot = uppot(children,:);
            childrenpot = reshape(childrenpot',[],1);
            
%             RHS = s(childrensurf,downsurf,epsilon,mu)*childrenpot;
            RHS = blockcomputation(childrensurf,downsurf,epsilon,mu,childrenpot);
            
        end
        
        pot = s(upsurf,downsurf,epsilon,mu) \ RHS;

        uppottemp(i,:) = pot';
    end
    
    uppot(nodes,:) = uppottemp;
    
end
toc

tic
for level = 1:1:levels
    nodes = find(tree.nodeLevel==level);

    downpottemp = zeros(size(nodes,2),coronapoints*3);

    for i = 1:size(nodes,2)
        node = nodes(i);
        
        upsurf = genupsurf(tree,node,coronares,coronarings);
        upsurf = reshape(upsurf.',[],1);

        downsurf = gendownsurf(tree,node,coronares,coronarings);
        downsurf = reshape(downsurf.',[],1);
        
        v = tree.interactions{node,2};
        x = tree.interactions{node,4};
        
        RHS = zeros(coronapoints*3,1);
        
        if ~isempty(v)
            
            vsurf = zeros(coronapoints*size(v,2),3);
            for vi = 1:size(v,2)
                vnode = v(vi);
                vsurf(((vi-1)*coronapoints)+1:vi*coronapoints,:) = ...
                    genupsurf(tree,vnode,coronares,coronarings);
            end
            vsurf = reshape(vsurf',[],1);
            
            vpot = reshape(uppot(v,:)',[],1);
           
%             RHS  = RHS + s(vsurf,upsurf,epsilon,mu)*vpot;
            RHS = RHS + blockcomputation(vsurf,upsurf,epsilon,mu,vpot);
            
        end
        
        if ~isempty(x)
            
            xsurf = zeros(tree.nodeCapacity*size(x,2),3);
            xpot = zeros(tree.nodeCapacity*size(x,2),3);
            for xi = 1:size(x,2)
                xnode = x(xi);

                xislice = tree.pointIndex == xnode;
                
                xipoints = tree.points(xislice,:);
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
            
%             RHS  = RHS + s(xsurf,upsurf,epsilon,mu)*xpot;
            RHS = RHS + blockcomputation(xsurf,upsurf,epsilon,mu,xpot);
            
        end
        
        parentnode = tree.nodeParents(node);
        
        parentsurf = gendownsurf(tree,parentnode,coronares,coronarings);
        parentsurf = reshape(parentsurf',[],1);
        
        parentpot = downpot(parentnode,:)';
        
%         RHS = RHS + s(parentsurf,upsurf,epsilon,mu) * parentpot;
        RHS = RHS + blockcomputation(parentsurf,upsurf,epsilon,mu,parentpot);
        
        pot = s(downsurf,upsurf,epsilon,mu) \ RHS;

        downpottemp(i,:) = pot';
    end
    
    downpot(nodes,:) = downpottemp;
    
end
toc

leaves = find(isleaf(tree,1:tree.nodeCount));

vel = zeros(size(potentials));
velpar = cell(size(leaves,1),1);

tic
for i = 1:size(leaves,1)
    
    node = leaves(i);
    
    downsurf = gendownsurf(tree,node,coronares,coronarings);
    downsurf = reshape(downsurf',[],1);
    
    nodepot = downpot(node,:)';
    
    leafslice = tree.pointIndex==node;
    leafpoints = tree.points(leafslice,:);
    leafpoints = reshape(leafpoints',[],1);
    
    veltemp = s(downsurf,leafpoints,epsilon,mu)*nodepot;
%     veltemp = blockcomputation(downsurf,leafpoints,epsilon,mu,nodepot);
    
    u = tree.interactions{node,1};
    w = tree.interactions{node,3};
    
    if ~isempty(u)

        usurf = zeros(tree.nodeCapacity*size(u,2),3);
        upot = zeros(tree.nodeCapacity*size(u,2),3);
        for ui = 1:size(u,2)
            unode = u(ui);

            uislice = tree.pointIndex == unode;

            uipoints = tree.points(uislice,:);
            uipointpot = potentials(uislice,:);

            points = size(uipoints,1);

            usurf(((ui-1)*tree.nodeCapacity)+1:...
                  ((ui-1)*tree.nodeCapacity)+points,:) = ...
                  uipoints;
            upot(((ui-1)*tree.nodeCapacity)+1:...
                 ((ui-1)*tree.nodeCapacity)+points,:) = ...
                 uipointpot;

        end
        usurf(~any(usurf,2),:) = [];
        upot(~any(upot,2),:) = [];

        usurf = reshape(usurf',[],1);
        upot = reshape(upot',[],1);

%         veltemp  = veltemp + s(usurf,leafpoints,epsilon,mu)*upot;
        veltemp = veltemp + blockcomputation(usurf,leafpoints,epsilon,mu,upot);

    end
    
    if ~isempty(w)
            
        wsurf = zeros(coronapoints*size(w,2),3);
        for wi = 1:size(w,2)
            wnode = w(wi);
            wsurf(((wi-1)*coronapoints)+1:wi*coronapoints,:) = ...
                genupsurf(tree,wnode,coronares,coronarings);
        end
        wsurf = reshape(wsurf',[],1);

        wpot = reshape(uppot(w,:)',[],1);

%         veltemp  = veltemp + s(wsurf,leafpoints,epsilon,mu)*wpot;
        veltemp = veltemp + blockcomputation(wsurf,leafpoints,epsilon,mu,wpot);

    end
    
    velpar{i} = reshape(veltemp,3,[])';
    
end

for i = 1:size(leaves,1)
    
    node = leaves(i);
    
    leafslice = tree.pointIndex==node;

    vel(leafslice,:) = velpar{i};
    
end
toc

diff = nystromvel-vel;
norm(nystromvel-vel,'fro')/norm(nystromvel,'fro')