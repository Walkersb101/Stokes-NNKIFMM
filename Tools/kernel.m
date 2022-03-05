function [S] = kernel(x,x0,param)
% Kernel  Compute the stokeslet matrix with x and x0
% Compute the a (3*M, 3*N) matrix based on the form of the regualraised
% stokeslet given in https://aip.scitation.org/doi/abs/10.1063/1.1830486
%Inputs:
%   x    : A (3*N,1) vector of stokelet points where we have the potential
%          formated as [x1 y1 z1 x2 y2 z2 ...]'
%   x0   : A (3*M,1) vector of sample points where we want to compute 
%          the velocity formated as [x1 y1 z1 x2 y2 z2 ...]'
%  param : Parameters for kernel
%
%Output:
%   S : (3*M, 3*N) stokeslet matrix 

epsilon = param(1);
mu = param(2);

% find number of partices in x and x0
N=length(x)/3;
M=length(x0)/3;


rx = (x0(1:3:end)-x(1:3:end)');
ry = (x0(2:3:end)-x(2:3:end)');
rz = (x0(3:3:end)-x(3:3:end)');
% generate r^2 for each particle
r2 = rx.^2 + ry.^2 + rz.^2;
% inverse ry with blob and expand to matrix size 
invs = kron((1./((r2+epsilon^2).^(3/2))).*(1/(8*pi*mu)),ones(3));
% compute first term in S
if isgpuarray(x)
    dyadic = zeros(3*M,3*N,'gpuArray');
else
    dyadic = zeros(3*M,3*N);
end
dyadic(1:3:end,1:3:end) = rx.*rx;
dyadic(1:3:end,2:3:end) = rx.*ry;
dyadic(1:3:end,3:3:end) = rx.*rz;
dyadic(2:3:end,1:3:end) = ry.*rx;
dyadic(2:3:end,2:3:end) = ry.*ry;
dyadic(2:3:end,3:3:end) = ry.*rz;
dyadic(3:3:end,1:3:end) = rz.*rx;
dyadic(3:3:end,2:3:end) = rz.*ry;
dyadic(3:3:end,3:3:end) = rz.*rz;

% compute second term in S
isotropic = kron(r2+2*epsilon^2,eye(3));
% Add terms and multiply by invs
S = (dyadic + isotropic).*invs;
clearvars isotropic dyadic
end