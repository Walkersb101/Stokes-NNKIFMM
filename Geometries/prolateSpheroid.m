function [points] = prolateSpheroid(a,c,h)
% prolateSpheroid  Discretize a prolate Spheroid
% Discretize a prolate spheroid with major axis c and minor axis a, the
% discretiation keeps the spacing of points to be approximatly h
%Inputs:
%   a      : Scalar value denoting the radius of the minor axis
%   c      : Scalar value denoting the radius of the major axis
%   centre : a (1,3) array denoting the centre of the spheroid
%   h      : approximate point spacing
%
%Output:
%   points   : (M, 3) array of point on the sphere

points = [0 0 c; 0 0 -c];

vN = ceil(pi*(3*(a+c)-sqrt((3*a+c)*(a+3*c)))/(2*h));
v = linspace(0,pi,vN);

for i = 2:vN-1
    z = (c*cos(v(i)))';
    uN = ceil(2*pi*a*sqrt(1-(z/c)^2)/h);
    u = linspace(0,2*pi,uN+1);
    u = u(2:end);
    x = (a*cos(u)*sin(v(i)))';
    y = (a*sin(u)*sin(v(i)))';
    z = z*ones(uN,1);
    points = vertcat(points, [x,y,z]);
end
end
