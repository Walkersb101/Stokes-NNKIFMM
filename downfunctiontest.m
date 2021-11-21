clear all

points = [Sphere(1,[-10 -10 -10],30);
          Sphere(1,[10 -10 -10],30);
          Sphere(1,[-10 10 -10],30);
          Sphere(1,[10 10 -10],30);
          Sphere(1,[-10 -10 10],30);
          Sphere(1,[10 -10 10],30);
          Sphere(1,[-10 10 10],30);
          Sphere(1,[10 10 10],30)];    
      
     
pointpot = 10*rand(size(points))-5;

resolution = 6;
shells = 2;

tree = OcTree(points,pointpot,'nodeCapacity',400,'maxDepth',21);

points = [];

interactions = geninteractionlists(tree);

index = 404;

coords = tree.nodeCorners(index,:);
width = coords(4:6)- coords(1:3);
nearboundary = [coords(1:3) - width*1.1 coords(4:6) + width*1.1];

domaincoords = tree.nodeCorners(1,:);
domainwidth = domaincoords(4:6)- domaincoords(1:3);
maxdepth = max(tree.nodeLevel);
smallestnode = domainwidth.*(2^(-maxdepth));

farboundary = [nearboundary(1:3)-smallestnode, ...
               nearboundary(4:6)+smallestnode];

interres = resolution - 2*(shells-1);
spacing = smallestnode/(shells-1);

x = [farboundary(1):spacing(1):nearboundary(1)-spacing(1),...
     linspace(nearboundary(1),nearboundary(4),interres),...
     farboundary(4):-spacing(1):nearboundary(4)+spacing(1)];
y = [farboundary(2):spacing(2):nearboundary(2)-spacing(2),...
     linspace(nearboundary(2),nearboundary(5),interres),...
     farboundary(5):-spacing(2):nearboundary(5)+spacing(2)];
z = [farboundary(3):spacing(3):nearboundary(3)-spacing(3),...
     linspace(nearboundary(3),nearboundary(6),interres),...
     farboundary(6):-spacing(3):nearboundary(6)+spacing(3)];

[X,Y,Z] = meshgrid(x,y,z);

X(shells+1:end-shells,shells+1:end-shells,shells+1:end-shells)= NaN;
Y(shells+1:end-shells,shells+1:end-shells,shells+1:end-shells)= NaN;
Z(shells+1:end-shells,shells+1:end-shells,shells+1:end-shells)= NaN;

X=reshape(X(~isnan(X)),[],1);
Y=reshape(Y(~isnan(Y)),[],1);
Z=reshape(Z(~isnan(Z)),[],1);


hold on
plot3(X, Y, Z,'o')

[U,V,W,X] = interactions{index,1:4};

plotcube(tree.nodeCorners(U,:),'g')
plotcube(tree.nodeCorners(V,:),'c')
plotcube(tree.nodeCorners(W,:),'k')
plotcube(tree.nodeCorners(X,:),'b')
plotcube(tree.nodeCorners(index,:),'r')