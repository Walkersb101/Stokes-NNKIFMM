function [points] = Rotate(points, alpha, beta, gamma, centre)
% Rotate  Rotate points about axis centred on centre
%Inputs:
%   alpha  : Angle about x axis
%   beta   : Angle about y axis
%   gamma  : Angle about z axis
%   centre : a (1,3) array denoting the centre rotation
%
%Output:
%   points   : (M, 3) array of point on the sphere

points = points - centre.*ones(size(points,1),1);

Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
rotation =  Rx*Ry*Rz;

points = (rotation * points')';

points = points + centre.*ones(size(points,1),1);
end

