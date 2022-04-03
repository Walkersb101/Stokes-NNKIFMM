function [points] = cylinderBoundary(r,height,centre,h)
%CYLINDERBOUNDARY Discretize a cylinder parallel to x axis
%Inputs:
%   r      : Radius of the cylinder
%   height : Height of the cylinder
%   centre : (1,3) array denoting centre of the cylinder
%   h      : Spacing between points
%
%Output:
%   points   : (N, 3) array of point on the plane

angleRes = ceil((2*pi*r)/h);
heightRes = ceil(height/h)+1;

theta = linspace(0,2*pi,angleRes+1).';
theta = repmat(theta(2:end),heightRes,1);
r = r*ones(size(theta));
z = kron(linspace(-height/2,height/2,heightRes).',ones(angleRes,1));

x = r .* cos(theta);
y = r .* sin(theta);

points = [x,y,z] + centre.*ones(angleRes*heightRes,1);
end

