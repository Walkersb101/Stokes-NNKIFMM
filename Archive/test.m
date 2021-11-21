clear all;
clf;

points = rand(400,3)*0.5;
points =  [points; 1 1 1];

pointpot = s(reshape(points.',[],1),reshape(points.',[],1),0.01,1) \ repmat([1;0;0],401,1);
pointpot = reshape(pointpot,3,[])';


tic
tree = OcTree(points,'binCapacity',200,'maxDepth',21);
toc

plot3(points(:,1), points(:,2), points(:,3),'.')
hold on;
axis equal;

tic
interactions = geninteractionlists(tree);
toc

upsurf = genupsurf(tree,10,6,2);
upsurf = reshape(upsurf.',[],1);

downsurf = gendownsurf(tree,10,6,2);
downsurf = reshape(downsurf.',[],1);

upsurf-downsurf


det(s(upsurf,downsurf,0.01,1))
