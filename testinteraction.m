clear all;
clf;

points = rand(10000,3);

tic
tree = OcTree(points,'binCapacity',100,'maxDepth',21);
toc

plot3(points(:,1), points(:,2), points(:,3),'.')
hold on;
axis equal;

tic
interactions = geninteractionlists(tree);
toc

randNodeIndex = randi([1 tree.binCount]);
[U,V,W,X] = interactions{randNodeIndex,1:4};



% Plot cubes
plotcube(tree.binCorners(U,:),'g')
plotcube(tree.binCorners(V,:),'c')
plotcube(tree.binCorners(W,:),'k')
plotcube(tree.binCorners(X,:),'b')
plotcube(tree.binCorners(randNodeIndex,:),'r')
hold off;
