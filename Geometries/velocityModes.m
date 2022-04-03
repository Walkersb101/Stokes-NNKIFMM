function [vel] = velocityModes(points,vel,mom)
%TRANSLATIONALMODE Compute velocity of points given translation modes
%
% Inputs:
%   points : (N,3) array of point positions
%   vel    : (1,3) array of rigid body velocity
%   mom    : (1,3) array of rigid body angular velocity
%
% Output:
%   vel : (N,3) array of point velocities


N = size(points,1);

x = points(:,1);
y = points(:,2);
z = points(:,3);

ux = vel(1)*ones(N,1);
uy = vel(2)*ones(N,1);
uz = vel(3)*ones(N,1);

point_vel = [ux; uy; uz] + [mom(2) * z - mom(3) * y; mom(3) * x - mom(1) * z; mom(1) * y - mom(2) * x];

vel = reshape(point_vel,[],3);
end

