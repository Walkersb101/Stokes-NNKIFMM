points = rand(100000,3);

tic
tree = OcTree(points,'binCapacity',10,'maxDepth',10);
toc



% hold on;
% axis equal;
% plot3(points(:,1), points(:,2), points(:,3),'.')
% 
% corners = tree.binCorners;
% 
% for i = 1:size(tree.binCorners,1)
%     Plotcube(corners(i,:))
% end