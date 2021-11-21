function [RHS] = blockcomputation(x,x0,epsilon,mu,potentials)
%BLOCKCOMPUTATION Summary of this function goes here
%   Detailed explanation goes here

x = gpuArray(x);
x0 = gpuArray(x0);
potentials = gpuArray(potentials);

N = size(x0,1)/3;

Max_array_size = (0.1 * 8000000000)/64;
Row_block = floor(Max_array_size/(N*3));
blocks = [1:3*Row_block:3*N (3*N)+1];

RHS = zeros(3*N,1);
for i = 1:length(blocks)-1
    RHS(blocks(i):blocks(i+1)-1,:) = s(x,x0(blocks(i):blocks(i+1)-1,:),epsilon,mu)*potentials;
end
RHS = gather(RHS);
end

