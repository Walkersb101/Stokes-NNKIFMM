function [points] = Sphere(r,centre, N)
% s  Discretize a sphere
% Discretize a sphere based on a cubic discretization. The sphere has
% radius r and centre centre. The cubic discretization has N points in each
% axis
%Inputs:
%   r      : Scalar value denoting the radius of the sphere
%   centre : a (1,3) array denoting the centre of the sphere
%   N      : Number of points in each axis of the discretization
%
%Output:
%   s   : (M, 3) array of point on the sphere


% generate square grid
pos = linspace(-1*r,r,N);
[X,Y,Z] = meshgrid(pos,pos,pos);

% remove centre of grid
X(2:end-1,2:end-1,2:end-1) = NaN;
Y(2:end-1,2:end-1,2:end-1) = NaN;
Z(2:end-1,2:end-1,2:end-1) = NaN;


X=(X(~isnan(X)));
Y=(Y(~isnan(Y)));
Z=(Z(~isnan(Z)));

% convert to polar
azimuth = atan2(Y,X);
elevation = atan2(Z,sqrt(X.^2 + Y.^2));

% convert to cartesian with new radius
x = r .* cos(elevation) .* cos(azimuth);
y = r .* cos(elevation) .* sin(azimuth);
z = r .* sin(elevation);

% move to centre
x = x + centre(1);
y = y + centre(2);
z = z + centre(3);

points = [x y z];

end

