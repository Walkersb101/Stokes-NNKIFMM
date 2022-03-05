function [x1, x2, x3] = extractComponents(x)
%extractComponents extracts componets of interleaved vector
%Takes a vector of interleaved componets and returns the 3 vectors which
%make it up
% Inputs:
%   x : A (3*N,1) vector of interleaves components
%
% Outputs:
%   x1 : A (N,1) or (1,N) vector
%   x2 : A (N,1) or (1,N) vector
%   x2 : A (N,1) or (1,N) vector

N = numel(x);
x1 = x(1:3:N);
x2 = x(2:3:N);
x3 = x(3:3:N);
end

