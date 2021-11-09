esp = 0.001;
mu = 1;

particlenumber = 10000;

points = rand(particlenumber,3);
u = rand(particlenumber,3);

tic
tree = OcTree(points,'binCapacity',200,'maxDepth',21);
toc

points = [];

interactions = geninteractionlists(tree);

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
    
    coords = tree.binCorners(nonleafi,:);
    upsurf = gencorona(coords,6,2);
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