function [B,BT,A1,A2] = genMatrixComponets(swimmers)
%genMatrixComponets Create matrices B, BT and A for computation in mobility
% problem, Reduced overall complexity as these matrices are needed multiple
% times and remain constant for each time step.
%
% Inputs:
%   swimmers : A Swimmer class construction containing swimmers and
%               boundary
%             
% Outputs:
%   B  : Lower matrix needed to compute total force and torque
%   BT : Upper matrix need to compute the velocity of individual
%        force points from the total velocity and angular velocity
%   A  : inverse of B*BT

Nearest = swimmers.Tree.arguments.nearest;
[~,~,NN,points,finePoints] = swimmers.getSwimmers();

Nsw = size(points,1);
points = reshape(points.',[],1);

bndPoints = swimmers.getBnd();
NBnd = size(bndPoints,1);

N = Nsw+NBnd;

noSwimmers = swimmers.swimmerNo();
swSizes = swimmers.getSwimmerSizes();
swIndex = zeros(1,noSwimmers+1);swIndex(1) = 1;
for i = 2:noSwimmers+1
    swIndex(i) = swIndex(i-1) + swSizes(i-1);
end

au = sparse(N,noSwimmers);
for n = 1 : noSwimmers
    au(swIndex(n):swIndex(n+1)-1,n) = -ones(swIndex(n+1)-swIndex(n),1);
end
AU = kron(speye(3),au);

Ze  = sparse(N, noSwimmers);
x1m = sparse(N, noSwimmers);
x2m = sparse(N, noSwimmers);
x3m = sparse(N, noSwimmers);
for n = 1 : noSwimmers
    swimPoints = swimmers.getIndSwimmer(n);
    x1m(swIndex(n):swIndex(n+1)-1,n) = swimPoints(:,1);
    x2m(swIndex(n):swIndex(n+1)-1,n) = swimPoints(:,2);
    x3m(swIndex(n):swIndex(n+1)-1,n) = swimPoints(:,3);
end
AOm = [ Ze -x3m x2m; x3m Ze -x1m; -x2m x1m Ze];


if Nearest
    Q = size(finePoints,1);
    finePoints = reshape(finePoints,[],1);
    NN = kron(speye(3),NN);
    
    af = sparse(noSwimmers, N);
    for n= 1 : noSwimmers
        af(n,swIndex(n):swIndex(n+1)-1) = ...
            sum(NN(1:Q,swIndex(n):swIndex(n+1)-1),1);
    end
    AF = kron(speye(3),af);
    
    Ze  = sparse(noSwimmers,N);
    x1m = sparse(noSwimmers,N);
    x2m = sparse(noSwimmers,N);
    x3m = sparse(noSwimmers,N);
    pointMap = finePoints'*NN;
    x1 = pointMap(1:Nsw);x2 = pointMap(Nsw+1:2*Nsw);x3 = pointMap(2*Nsw+1:end);
    for n=1:noSwimmers
        x1m(n,swIndex(n):swIndex(n+1)-1)=x1(swIndex(n):swIndex(n+1)-1);
        x2m(n,swIndex(n):swIndex(n+1)-1)=x2(swIndex(n):swIndex(n+1)-1);
        x3m(n,swIndex(n):swIndex(n+1)-1)=x3(swIndex(n):swIndex(n+1)-1);
    end
    AM = [Ze -x3m x2m; x3m Ze -x1m; -x2m x1m Ze];
else
    af = sparse(noSwimmers, N);
    for n= 1 : noSwimmers
        af(n,swIndex(n):swIndex(n+1)-1) = ones(swIndex(n+1)-swIndex(n),1);
    end
    AF = kron(speye(3),af);

    Ze  = sparse(noSwimmers,N);
    x1m = sparse(noSwimmers,N);
    x2m = sparse(noSwimmers,N);
    x3m = sparse(noSwimmers,N);
    [x1,x2,x3]=extractComponents(points);
    for n=1:noSwimmers
        x1m(n,swIndex(n):swIndex(n+1)-1)=x1(swIndex(n):swIndex(n+1)-1);
        x2m(n,swIndex(n):swIndex(n+1)-1)=x2(swIndex(n):swIndex(n+1)-1);
        x3m(n,swIndex(n):swIndex(n+1)-1)=x3(swIndex(n):swIndex(n+1)-1);
    end
    AM = [Ze -x3m x2m; x3m Ze -x1m; -x2m x1m Ze];
end

B = [AF;AM];
BT = [AU AOm];
A1 = B*B.';
A2 = BT.'*BT;
end

