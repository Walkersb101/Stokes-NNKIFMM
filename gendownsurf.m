function [surf] = gendownsurf(tree, index, resolution, shells)
%GENUPSURFACE Summary of this function goes here
%   Detailed explanation goes here


    

    coords = tree.binCorners(index,:);
    width = coords(4:6)- coords(1:3);
    nearboundary = [coords(1:3) - width coords(4:6) + width];
    
    interiorres = resolution-2*(shells-1);
    spacing = 3*width/interiorres;
    coronadepth = (shells-1)*spacing;
    
    outerboundary = [nearboundary(1:3)-coronadepth, ...
                     nearboundary(4:6)+coronadepth];
    
    surf = gencorona(outerboundary,resolution,shells);
end

