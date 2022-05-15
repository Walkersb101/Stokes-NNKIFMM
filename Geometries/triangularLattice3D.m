function [X,Y,Z] = triangularLattice3D(h_dist,centre,limits)
%TRAINGULARLATTICE Form triangualar lattice
%   Compute triangular lattice on YZ Plane. Horizontal discretiasion is
%   h_dist and will tile until to limits. Lattice is then translated to
%   centre
%
% Input:
%   h_dist : horizontal discretisation distance
%   centre : [X,Y,Z] centre
%   limits : size of tiling [xLimit,yLimit,zLimit]
%
% Output:
%   X : (N,1) array of X coordinates
%   Y : (N,1) array of Y coordinates

v_dist = sqrt(h_dist^2-(h_dist/2)^2);

x_lim = limits(1);
y_lim = limits(2);
z_lim = limits(3);

trigrid = [];
y_current = 0;
z_current = 0;
countz = 0;
while z_current < z_lim
    if mod(countz,2) == 0
        countx = 0;
        y_current = 0;
    else
        countx = 1;
        y_current = h_dist/2;
    end
    while y_current < y_lim
        if mod(countx,2) == 0
            xx = (0:h_dist:x_lim)';
            yy = ones(length(xx), 1)*y_current;
            zz = ones(length(xx), 1)*z_current;
            displacementx = 1;
        else
            xx = (h_dist/2:h_dist:x_lim)';
            yy = ones(length(xx), 1)*y_current;
            zz = ones(length(xx), 1)*z_current;
            displacementx = 0;
        end
        trigrid = [trigrid; [xx,yy,zz]];
        y_current = y_current + v_dist;
        countx = countx + 1;
    end
    z_current = z_current + v_dist;
    countz = countz + 1;
end
X = trigrid(:,1)-(x_lim/2)+centre(1); Y = trigrid(:,2)-(y_lim/2)+centre(2); Z = trigrid(:,3)-(z_lim/2)+centre(3);
end

