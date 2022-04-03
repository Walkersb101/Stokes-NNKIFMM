function [points] = wall(centre,scale,axis,h)
% wall  Discretize a plane parallel to an axis
% Discretize a plane based on a linear discretization.
%Inputs:
%   centre : A (1,3) array denoting the centre of the plane
%   scale  : A (1,2) array denoting the size in each axis
%   axis   : 1 for parallel to X axis, 2 for parallel to Y axis and
%                2 for parallel to Z axis
%   h      : Horizontal and vertical spacing between points
%
%Output:
%   points   : (N, 3) array of point on the plane

resolution = ceil(scale/h);
noPoints = prod(resolution);
if axis == 1
    y = linspace(centre(2)-scale(1),centre(2)+scale(1),resolution(1));
    z = linspace(centre(3)-scale(2),centre(3)+scale(2),resolution(2));
    [Y,Z] = meshgrid(y,z);
    points = [centre(1)*ones(noPoints,1),Y(:),Z(:)];
elseif axis == 2
    x = linspace(centre(1)-scale(1),centre(1)+scale(1),resolution(1));
    z = linspace(centre(3)-scale(2),centre(3)+scale(2),resolution(2));
    [X,Z] = meshgrid(x,z);
    points = [X(:),centre(2)*ones(noPoints,1),Z(:)];
else 
    x = linspace(centre(1)-scale(1),centre(1)+scale(1),resolution(1));
    y = linspace(centre(2)-scale(2),centre(2)+scale(2),resolution(2));
    [X,Y] = meshgrid(x,y);
    points = [X(:),Y(:),centre(3)*ones(noPoints,1)];
end
end

