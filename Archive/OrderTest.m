clear all;

points = rand(10000,3);

tic
tree = OcTree(points,'binCapacity',100,'maxDepth',21);
toc

tree.binCount
postordertree(tree,1)

preordertree(tree,1)