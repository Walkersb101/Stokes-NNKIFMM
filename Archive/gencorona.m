function [points] = gencorona(coords, resolution, shells)
%GENERATECORONA Generate square corona with set resolution and rings
%   Generate a square corna with discritisation of
%   (resolution,resolution,resolution) and corners given by coords. Shells
%   determines the number of outer layers in the corona (1 gives a surface)
%Inputs:
%   coords     : a (1,6) matrix where [x1 y1 z1 x2 y2 z2] are the corners 
%                of the corona. [x1 y1 z1] denotes the lower left corner 
%                and [x2 y2 z2] denotes the upper right corner.
%   resolution : Number of discritisation points in each axis
%   rings      : Number of outer shells
%
%Output:
%   points     : (n,3) array of points on the outer shell

x = linspace(coords(1),coords(4),resolution);
y = linspace(coords(2),coords(5),resolution);
z = linspace(coords(3),coords(6),resolution);

[X,Y,Z] = meshgrid(x,y,z);

X(shells+1:end-shells,shells+1:end-shells,shells+1:end-shells)= NaN;
Y(shells+1:end-shells,shells+1:end-shells,shells+1:end-shells)= NaN;
Z(shells+1:end-shells,shells+1:end-shells,shells+1:end-shells)= NaN;

X=reshape(X(~isnan(X)),[],1);
Y=reshape(Y(~isnan(Y)),[],1);
Z=reshape(Z(~isnan(Z)),[],1);

points = [X Y Z];

end

