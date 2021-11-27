function [s] = kernel(x,x0,param)
% s  Compute the stokeslet matrix with x and x0
% Compute the a (3*M, 3*N) matrix based on the form of the regualraised
% stokeslet given in https://aip.scitation.org/doi/abs/10.1063/1.1830486
%Inputs:
%   x    : A (3*N,1) vector of stokelet points where we have the force
%          formated as [x1 y1 z1 x2 y2 z2 ...]'
%   x0   : A (3*M,1) vector of sample points where we want to compute 
%          the velocity formated as [x1 y1 z1 x2 y2 z2 ...]'
%  param : Parameters for kernel
%
%Output:
%   s : (3*M, 3*N) stokeslet matrix 

epsilon = param(1);
mu = param(2);

% find number of partices in x and x0
N=length(x)/3;
M=length(x0)/3;

% generate r^2 for each particle
r2 = (x0(1:3:end)-x(1:3:end)').^2+(x0(2:3:end)-x(2:3:end)').^2+...
    (x0(3:3:end)-x(3:3:end)').^2;
% inverse r2 with blob and expand to matrix size 
invs = kron(1./((r2+epsilon^2).^(3/2)),ones(3)).*(1/(8*pi*mu)); 

% compute first term in S
dyadic = (repmat(kron(reshape(x,3,[]),ones(1,3)),M,1)-x0.*ones(1,3*N)).*...
    (ones(3*M,1).*x'-kron(repmat(reshape(x0,3,[])',1,N),ones(3,1)));

% compute second term in S
isotropic = kron(r2+2*epsilon^2,eye(3));
% Add terms and multiply by invs
s = (dyadic + isotropic).*invs; 
end