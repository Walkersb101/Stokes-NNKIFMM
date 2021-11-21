clear all;

epsilon = 0.001;
mu = 1;

% points = rand(400,3)*0.5;
% points =  [points; 1 1 1];
% 
% pointpot = rand(401,3);

points = rand(2000,3);
pointpot = rand(2000,3);

tic

N = size(points,1);
pointsreshape = reshape(points.',[],1);
pointpotreshape = reshape(pointpot.',[],1);

Max_array_size = (0.2 * 8000000000)/64;
Row_block = floor(Max_array_size/(N*3));
blocks = [1:3*Row_block:3*N (3*N)+1];

nystromvel = zeros(3*N,1);
for i = 1:length(blocks)-1
    nystromvel(blocks(i):blocks(i+1)-1,:) = s(pointsreshape,pointsreshape(blocks(i):blocks(i+1)-1,:),epsilon,mu)*pointpotreshape;
end
nystromvel = reshape(nystromvel,3,[])';
toc


tic
tree = OcTree(points,pointpot'binCapacity',200,'maxDepth',21);

points = [];

interactions = geninteractionlists(tree);
toc

postorder = postordertree(tree,1);
postorderslice = isleaf(tree,postorder);
postorderleaves = postorder(postorderslice);
postordernonleaves = postorder(~postorderslice);

uppot = cell(tree.binCount,1);

tic
for i = 1:size(postorderleaves,2)
    
    leafi = postorderleaves(i);
    
    upsurf = genupsurf(tree,leafi,6,2);
    upsurf = reshape(upsurf.',[],1);
    
    leafislice = find(tree.pointIndex == leafi);
    leafipoints = tree.points(leafislice,:);
    leafipoints = reshape(leafipoints.',[],1);
    
    leafipot = pointpot(leafislice,:);
    leafipot = reshape(leafipot.',[],1);
    
    vel = s(leafipoints,upsurf,epsilon,mu) * leafipot;
    
    pot = s(upsurf,upsurf,epsilon,mu) \ vel;
    uppot{leafi} = pot;
end

for i = 1:size(postordernonleaves,2)
    
    nonleafi = postordernonleaves(i);
    
    upsurf = genupsurf(tree,nonleafi,6,2);
    upsurf = reshape(upsurf.',[],1);
    
    children = tree.binChildren(nonleafi,:);
    childrenuppot = cat(1,uppot{children});
    
    childrenupsurf = [];
    for j = 1:8
        childrenupsurf = [childrenupsurf ; ...
                          reshape(genupsurf(tree,children(j),6,2).',[],1)];
    end
    
    RHS = s(childrenupsurf,upsurf,epsilon,mu)*childrenuppot;
    
    pot = s(upsurf,upsurf,epsilon,mu) \ RHS;
    
    uppot{nonleafi} = pot;
end
toc

downpot = cell(tree.binCount,1);

% upsurf = genupsurf(tree,1,6,2);
% upsurf = reshape(upsurf.',[],1);
% 
% pot = uppot{1};
% 
% downsurf = gendownsurf(tree,1,6,2);
% downsurf = reshape(downsurf.',[],1);
% 
% RHS = s(upsurf,downsurf,epsilon,mu)*pot;
% 
% downpot{1} = s(downsurf,downsurf,epsilon,mu) \ RHS;

preorder = preordertree(tree,1);
preorderleaves = preorder(isleaf(tree,preorder));

tic
for i = 2:size(preorder,2)
    
    nonrooti = preorder(i);
    
    coords = tree.binCorners(nonrooti,:);
    downsurf = gendownsurf(tree,nonrooti,6,2);
    downsurf = reshape(downsurf.',[],1);
    
    RHS = zeros(size(downsurf));
    
    v = interactions{nonrooti,2};
    x = interactions{nonrooti,4};
    
    if ~isempty(v)
        parfor vi = 1:size(v,2)
            vnode = v(vi);
            surf = genupsurf(tree,vnode,6,2);
            surf = reshape(surf.',[],1);
            pot = uppot{vi};
            RHS = RHS + s(surf,downsurf,epsilon,mu) * pot;
        end
    end
    
    if ~isempty(x)
        for xi = 1:size(x,2)
            xnode = x(xi);
            
            xislice = find(tree.pointIndex == xnode);
            xipoints = tree.points(xislice,:);
            xipoints = reshape(xipoints.',[],1);
            
            xipointpot = pointpot(xislice,:);
            xipointpot = reshape(xipointpot.',[],1);
            
            RHS = RHS + s(xipoints,downsurf,epsilon,mu) * xipointpot;
            
        end
    end
    
    pnode = tree.binParents(nonrooti);
    if pnode ~=1
        psurf = gendownsurf(tree,pnode,6,2);
        psurf = reshape(psurf.',[],1);
        ppot = downpot{pnode};
    
        RHS = RHS + s(psurf,downsurf,epsilon,mu) * ppot;
    end
    

    pot = s(downsurf,downsurf,epsilon,mu) \ RHS;
    
    downpot{nonrooti} = pot;
end

vel = zeros(size(tree.points));

for i = 1:size(preorderleaves,2)
    
    leafi = preorderleaves(i);
    
    leafislice = find(tree.pointIndex == leafi);
    leafipoints = tree.points(leafislice,:);
    leafipoints = reshape(leafipoints.',[],1);
    
    downsurf = gendownsurf(tree,leafi,6,2);
    downsurf = reshape(downsurf.',[],1);
    
    pot = downpot{leafi};
    
    pointvel = zeros(size(leafipoints));
    pointvel = s(downsurf,leafipoints,epsilon,mu) * pot;

    
    u = interactions{leafi,1};
    w = interactions{leafi,3};
   
    if ~isempty(u)
        for ui = 1:size(u,2)
            unode = u(ui);
            
            uislice = find(tree.pointIndex == unode);
            uipoints = tree.points(uislice,:);
            uipoints = reshape(uipoints.',[],1);
            
            uipointpot = pointpot(uislice,:);
            uipointpot = reshape(uipointpot.',[],1);
            
            pointvel = pointvel + s(uipoints,leafipoints,epsilon,mu) * uipointpot;
        end
    end
    
    if ~isempty(w)
        for wi = 1:size(w,2)
            wnode = w(wi);
            
            upsurf = genupsurf(tree,wnode,6,2);
            upsurf = reshape(upsurf.',[],1);
            
            pot = uppot{wnode};
            
            pointvel = pointvel + s(upsurf,leafipoints,epsilon,mu) * pot;
        end
    end
    
    vel(leafislice,:) = vel(leafislice,:) + reshape(pointvel,3,[])'; 
end
toc


diff = nystromvel-vel;
norm(nystromvel-vel,'fro')/norm(nystromvel,'fro')