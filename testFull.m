clear all;

esp = 0.001;
mu = 1;

particlenumber = 10000;

points = rand(particlenumber,3);
u = rand(particlenumber,3);

tic
tree = OcTree(points,'binCapacity',200,'maxDepth',21);

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
    
    leafivel = u(leafislice,:);
    leafivel = reshape(leafivel.',[],1);
    
    A = s(upsurf,leafipoints,esp,mu);
    
    pot = A \ leafivel;
    uppot{leafi} = pot;
end

for i = 1:size(postordernonleaves,2)
    
    nonleafi = postordernonleaves(i);
    
    upsurf = genupsurf(tree,leafi,6,2);
    upsurf = reshape(upsurf.',[],1);
    
    children = tree.binChildren(nonleafi,:);
    childrenuppot = cell2mat(uppot(children));
    
    childrenupsurf = [];
    for j = 1:8
        childrenupsurf = [childrenupsurf ; ...
                          reshape(genupsurf(tree,nonleafi,6,2).',[],1)];
    end
    
    RHS = s(childrenupsurf,upsurf,esp,mu)*childrenuppot;
    
    A = s(upsurf,upsurf,esp,mu);
    
    pot = A \ RHS;
    
    uppot{nonleafi} = pot;
end
toc

preorder = preordertree(tree,1);
preorderleaves = preorder(isleaf(tree,preorder));

downpot = cell(tree.binCount,1);

upsurf = genupsurf(tree,1,6,2);
upsurf = reshape(upsurf.',[],1);

pot = uppot{1};

downsurf = gendownsurf(tree,1,6,2);
downsurf = reshape(downsurf.',[],1);

RHS = s(upsurf,downsurf,esp,mu)*pot;

downpot{1} = s(downsurf,downsurf,esp,mu) \ RHS;

tic
for i = 2:size(preorder,2)
    
    nonrooti = preorder(i);
    
    coords = tree.binCorners(nonrooti,:);
    downsurf = gendownsurf(tree,nonrooti,6,2);
    downsurf = reshape(downsurf.',[],1);
    
    RHS = zeros(size(downsurf,2),1);
    
    v = interactions{nonrooti,2};
    x = interactions{nonrooti,4};
    
    if ~isempty(v)
        parfor vi = 1:size(v,2)
            vnode = v(vi);
            surf = genupsurf(tree,vnode,6,2);
            surf = reshape(surf.',[],1);
            pot = uppot{vi};
            RHS = RHS + s(surf,downsurf,esp,mu) * pot;
        end
    end
    
    if ~isempty(x)
        parfor xi = 1:size(x,2)
            xnode = x(xi);
            surf = genupsurf(tree,xnode,6,2);
            surf = reshape(surf.',[],1);
            pot = uppot{xi};
            RHS = RHS + s(surf,downsurf,esp,mu) * pot;
        end
    end
    
    pnode = tree.binParents(nonrooti);
    surf = genupsurf(tree,pnode,6,2);
    surf = reshape(surf.',[],1);
    pot = downpot{pnode};
    
    RHS = RHS + s(surf,downsurf,esp,mu) * pot;
    
    A = s(downsurf,downsurf,esp,mu);
    
    pot = A \ RHS;
    
    downpot{nonrooti} = pot;
end




toc

