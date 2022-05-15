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

[x1,x2] = triangularLattice(h,[0 0],scale);
noPoints = length(x1);
if axis == 1
    points = [centre(1)*ones(noPoints,1),x1-centre(2),x2-centre(3)];
elseif axis == 2
    points = [x1-centre(1),centre(2)*ones(noPoints,1),x2-centre(3)];
else 
    points = [x1-centre(1),x2-centre(2),centre(3)*ones(noPoints,1)];
end
end

