function [B,BT,A] = genMatrixComponets(swimmers)
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
[points,finePoints,NN] = swimmers.getSwimmers();

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

au = zeros(N,noSwimmers);
for n = 1 : noSwimmers
    au(swIndex(n):swIndex(n+1)-1,n) = -ones(swIndex(n+1)-swIndex(n),1);
end
AU = kron(eye(3),au);

Ze  = zeros(N, noSwimmers);
x1m = zeros(N, noSwimmers);
x2m = zeros(N, noSwimmers);
x3m = zeros(N, noSwimmers);
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
    
    af = zeros(noSwimmers, N);
    for n= 1 : noSwimmers
        af(n,swIndex(n):swIndex(n+1)-1) = ...
            sum(NN(1:Q,swIndex(n):swIndex(n+1)-1),1);
    end
    AF = kron(eye(3),af);
    
    Ze  = zeros(noSwimmers,N);
    x1m = zeros(noSwimmers,N);
    x2m = zeros(noSwimmers,N);
    x3m = zeros(noSwimmers,N);
    pointMap = finePoints'*NN;
    x1 = pointMap(1:Nsw);x2 = pointMap(Nsw+1:2*Nsw);x3 = pointMap(2*Nsw+1:end);
    for n=1:noSwimmers
        x1m(n,swIndex(n):swIndex(n+1)-1)=x1(swIndex(n):swIndex(n+1)-1);
        x2m(n,swIndex(n):swIndex(n+1)-1)=x2(swIndex(n):swIndex(n+1)-1);
        x3m(n,swIndex(n):swIndex(n+1)-1)=x3(swIndex(n):swIndex(n+1)-1);
    end
    AM = [Ze -x3m x2m; x3m Ze -x1m; -x2m x1m Ze];
else
    af = zeros(noSwimmers, N);
    for n= 1 : noSwimmers
        af(n,swIndex(n):swIndex(n+1)-1) = ones(swIndex(n+1)-swIndex(n),1);
    end
    AF = kron(eye(3),af);

    Ze  = zeros(noSwimmers,N);
    x1m = zeros(noSwimmers,N);
    x2m = zeros(noSwimmers,N);
    x3m = zeros(noSwimmers,N);
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
A = inv(B*BT);
end

