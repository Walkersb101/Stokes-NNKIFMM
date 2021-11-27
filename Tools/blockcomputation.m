function [RHS] = blockcomputation(x,x0,potentials,blockSize,kernalPar)
%BLOCKCOMPUTATION Compute matrix vector product of kernel with potential
%   Calcuated the matrix vector product of kernel with potential in block
%   of maximum size blockSize
% Inputs:
%   x          : A (3*N,1) vector of stokelet points where we have the 
%                force formated as [x1 y1 z1 x2 y2 z2 ...]'
%   x0         : A (3*M,1) vector of sample points where we want to compute 
%                the velocity formated as [x1 y1 z1 x2 y2 z2 ...]'
%   blockSize  : Maximum size of kernel matrix in GB
%   kernalPar  : Parameters for kernel 
%
% Outputs:
%   RHS : A (3*M,1) of velocities formated as 
%         [vx1 vy1 vZ1 vx2 vy2 vz2 ...]'

threeN = size(x0,1);

Max_array_size = (blockSize * 8000000000)/64;
Row_block = floor(Max_array_size/(threeN));
blocks = [1:3*Row_block:threeN threeN+1];

RHS = zeros(threeN,1);
for i = 1:length(blocks)-1
    RHS(blocks(i):blocks(i+1)-1,:) = ...
        kernel(x,x0(blocks(i):blocks(i+1)-1,:),kernalPar)*potentials;
end

end

