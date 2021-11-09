function [surf] = genupsurf(tree, index, resolution, shells)
%GENUPSURFACE Summary of this function goes here
%   Detailed explanation goes here

    coords = tree.binCorners(index,:); 
    surf = gencorona(coords,resolution,shells);
end

