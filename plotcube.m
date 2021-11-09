function plotcube(coords, varargin)
% plotcube  Plot cubes on the current plot
% plots cubes on the current plot based on the coords given
%Inputs:
%   coords   : a (n,6) matrix where [x1 y1 z1 x2 y2 z2] are the corners of
%              the cube. [x1 y1 z1] denotes the lower left corner 
%              and [x2 y2 z2] denotes the upper right corner.
%   varargin : variables for plot function
%

for i=1:size(coords,1)
    coord = coords(i,:);
    pts = cat(1, coord([...
                        1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                        1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                        nan(1,3), coord([4 2 3; 4 2 6]),...
                        nan(1,3), coord([4 5 3; 4 5 6]),...
                        nan(1,3), coord([1 5 3; 1 5 6]));
    plot3(pts(:,1),pts(:,2),pts(:,3), varargin{:});
end
