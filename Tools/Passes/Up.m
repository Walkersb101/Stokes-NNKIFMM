function [pot] = Up(nodepoints,upsurf,downsurf,nodepotentials,blockSize,kernelPar,GPU)
% Up  Compute the affect of points on the upwards equivlient surface
% Compute the effect of points on the downward equivilent potential
% before inverting to find the equilent potential on the upwards surface
%Inputs:
%   nodepoints     : A (3*N,1) vector of source points
%   upsurf         : A (3*Nq,1) vector defining the position of the upward 
%                    surface
%   downsurf       : A (3*Nq,1) vector defining the position of the 
%                    downward surface
%   nodepotentials : A (3*N,1) vector of source potentials
%   blockSize      : Block size of matrix calulations
%   kernelPar      : Paramaters for the kernel
%   GPU            : 1 for GPU compute, 0 otherwise
%
%Output:
%   pot : (1, 3*Nq) vector of potentials at the upward surface

if GPU
    nodepoints = gpuArray(nodepoints);
    nodepotentials = gpuArray(nodepotentials);
    upsurf = gpuArray(upsurf);
    downsurf = gpuArray(downsurf);
    
    pot = gather((kernel(upsurf,downsurf,kernelPar) \...
    blockcomputation(nodepoints,downsurf,nodepotentials,...
                     blockSize,kernelPar))');
                 
    clear nodepoints upsurf downsurf nodepotentials;
else
    pot = (kernel(upsurf,downsurf,kernelPar) \...
    blockcomputation(nodepoints,downsurf,nodepotentials,...
                     blockSize,kernelPar))';

end

