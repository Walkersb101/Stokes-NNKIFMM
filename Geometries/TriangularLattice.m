function [X,Y] = triangularLattice(h_dist,centre,limits)
%TRAINGULARLATTICE Form triangualar lattice
%   Compute triangular lattice on YZ Plane. Horizontal discretiasion is
%   h_dist and will tile until to limits. Lattice is then translated to
%   centre
%
% Input:
%   h_dist : horizontal discretisation distance
%   centre : [X,Y] centre
%   limits : size of tiling [xLimit,yLimit]
%
% Output:
%   X : (N,1) array of X coordinates
%   Y : (N,1) array of Y coordinates

v_dist = sqrt(h_dist^2-(h_dist/2)^2);

x_lim = limits(1);
y_lim = limits(2);

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
X = trigrid(:,1)-(x_lim/2)+centre(1); Y = trigrid(:,2)-(y_lim/2)+centre(2);
end

