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

% randNodeIndex = randi([1 tree.binCount]);
randNodeIndex = 9;
[U,V,W,X] = interactions{randNodeIndex,1:4};

% Plot cubes
plotcube(tree.binCorners(U,:),'g')
plotcube(tree.binCorners(V,:),'c')
plotcube(tree.binCorners(W,:),'k')
plotcube(tree.binCorners(X,:),'b')
plotcube(tree.binCorners(randNodeIndex,:),'r')

upsurf = genupsurf(tree,randNodeIndex,6,2);
plot3(upsurf(:,1),upsurf(:,2),upsurf(:,3),'.');

downsurf = gendownsurf(tree,randNodeIndex,6,2);
plot3(downsurf(:,1),downsurf(:,2),downsurf(:,3),'.');

hold off;
% 
% 
% leafi = 10;
% 
% upsurf = genupsurf(tree,leafi,6,2);
% upsurf = reshape(upsurf.',[],1);
% 
% downsurf = gendownsurf(tree,leafi,6,2);
% downsurf = reshape(downsurf.',[],1);
% 
% leafislice = find(tree.pointIndex == leafi);
% leafipoints = tree.points(leafislice,:);
% leafipoints = reshape(leafipoints.',[],1);
% 
% leafipot = pointpot(leafislice,:);
% leafipot = reshape(leafipot.',[],1);
% 
% vel = s(leafipoints,downsurf,0.01,1) * leafipot;
% 
% pot = s(upsurf,downsurf,0.01,1) \ vel;
% 
% 
% direct = s(reshape(leafipoints.',[],1),points(401,:)',0.01,1)*reshape(leafipot.',[],1)
% 
% indirect = s(upsurf,points(401,:)',0.01,1)*pot
% 
% norm(direct-indirect,'fro')/norm(direct,'fro')