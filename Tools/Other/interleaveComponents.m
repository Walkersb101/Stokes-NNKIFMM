function [x] = interleaveComponents(x1, x2, x3)
%interleaveComponents iterleaves components of x1, x2, x3
%Takes row or column vector, interleaves the componets and returns a new
%column vector
% Inputs:
%   x1 : A (N,1) or (1,N) vector
%   x2 : A (N,1) or (1,N) vector
%   x2 : A (N,1) or (1,N) vector
%
% Outputs:
%   x : A (3*N,1) vector of interleaves components

if isrow(x1)
   x = [x1;x2;x3];
   x = x(:)';
else
    x = [x1, x2, x3]';
    x = x(:);
end
end

