function [vel] = Downleaf(leafpoints,downsurf,downpot,usurf,upot,wsurf,wpot,blockSize,kernelPar,GPU)
% UpNonleaf  Compute the affect of nodes on the target points
%Inputs:
%   leafpoints : A vector of positions of target points in the node
%   downsurf   : A (3*Nq,1) vector defining the position of the 
%                downward surface
%   downpot    : A (3*Nq,1) vector defining the potential of the 
%                downward surface
%   usurf      : A vector of positions at u nodes upward surfaces
%   upot       : A vector of potentials at u nodes upward surfaces
%   wsurf      : A vector of positions at w nodes upward surfaces
%   wpot       : A vector of potentials at w nodes upward surfaces
%   kernelPar  : Paramaters for the kernel
%   GPU        : 1 for GPU compute, 0 otherwise
%
%Output:
%   vel : A vector of velocities at the target points

if GPU
    leafpoints = gpuArray(leafpoints);
    downsurf = gpuArray(downsurf);
    downpot = gpuArray(downpot);
    usurf = gpuArray(usurf);
    upot = gpuArray(upot);
    wsurf = gpuArray(wsurf);
    wpot = gpuArray(wpot);

    vel = blockcomputation(downsurf,leafpoints,downpot,blockSize,...
                           kernelPar)+...
          blockcomputation(usurf,leafpoints,upot,blockSize,kernelPar) +...
          blockcomputation(wsurf,leafpoints,wpot,blockSize,kernelPar);

    vel = gather(reshape(vel,3,[])');

    clear leafpoints downsurf nodepot usurf upot wsurf wpot
else
    vel = blockcomputation(downsurf,leafpoints,downpot,blockSize,...
                           kernelPar)+...
          blockcomputation(usurf,leafpoints,upot,blockSize,kernelPar) +...
          blockcomputation(wsurf,leafpoints,wpot,blockSize,kernelPar);
      
    vel = reshape(vel,3,[])';
end
end

