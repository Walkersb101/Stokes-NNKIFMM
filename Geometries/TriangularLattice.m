function [trigrid] = TriangularLattice(h_dist,limits)
%TRAINGULARLATTICE Summary of this function goes here
%   Detailed explanation goes here
%% Triangular grid information
v_dist = sqrt(h_dist^2-(h_dist/2)^2); % Vertical distance
%% Region size
x_lim = limits(1);
y_lim = limits(2);
%% Generate grid
trigrid = [];
y_current = 0;
xx = 0;
displacement = 0;
while y_current < y_lim
    if displacement == 0
        xx = [0:h_dist:x_lim]';
        yy = ones(length(xx), 1)*y_current;
        displacement = 1;
    else
        xx = [h_dist/2:h_dist:x_lim]';
        yy = ones(length(xx), 1)*y_current;
        displacement = 0;
    end
    trigrid = [trigrid; [xx,yy]];
    y_current = y_current + v_dist;
end
end

