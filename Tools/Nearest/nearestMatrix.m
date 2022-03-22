function [NN] = nearestMatrix(coarse,fine,blockSize,threads)


coarse = reshape(coarse.',[],1);
fine = reshape(fine.',[],1);

N = size(coarse,1)/3;
Q = size(fine,1)/3;

blockNodes=floor(blockSize*2^27/(9*Q));
blocks = [1:3*blockNodes:3*Q (3*Q)+1];

rows = [];
cols = [];
weight = [];
parfor (i = 1:length(blocks)-1, threads)
    d2 = sum((repmat(reshape(reshape(fine(blocks(i):blocks(i+1)-1)'...
              ,3,[])',[],1,3),1,N) - ...
              repmat(pagetranspose(reshape(reshape(coarse,3,[])',[],1,3))...
              ,(blocks(i+1)-blocks(i))/3,1)).^2,3);
    [A,~] = min(d2,[],2);
    [Ir,Ic] = find(d2 == A);
    values = ((d2 == A)./sum((d2 == A),2))';
    weight = [weight,values((d2 == A)')'];
    cols = [cols, Ic.'];
    rows = [rows, (blocks(i)-1)/3+Ir.'];
end
NN = sparse(rows,cols,weight,Q,N);
end

