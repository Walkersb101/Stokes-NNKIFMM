function [b] = MobilityMatrixMulti(f,swimmers,initalCon)

Nearest = swimmers.Tree.arguments.nearest;
points = swimmers.getSwimmerBnd();

noSwimmers = swimmers.swimmerNo();
swSizes = swimmers.getSwimmerSizes();
swIndex = zeros(1,noSwimmers+1);swIndex(1) = 1;
for i = 2:noSwimmers+1
    swIndex(i) = swIndex(i-1) + swSizes(i-1);
end

x0=initalCon(1:3);
b1=initalCon(4:6);
b2=initalCon(7:9);
b3=cross(b1,b2);
B=[b1(:) b2(:) b3(:)];

points = (B * points')';
points = points + x0'.*ones(size(points,1),1);

force = f(1:end-6*noSwimmers);
U = f(end-6*noSwimmers+1:end-3*noSwimmers);
Omega = f(end-3*noSwimmers+1:end);

N=numel(points)/3;

if Nearest
    [~, finepoints, NN] = swimmers.getSwimmerBnd();
    Q = numel(finepoints)/3;
    finepoints = (B * finepoints')';
    finepoints = finepoints + x0'.*ones(size(finepoints,1),1);
    finepoints = reshape(finepoints.',[],1);
end

b = zeros(3*N+6*noSwimmers,1);

swimmers.FMM.arguments.GMRES = 1;
b(1:end-6*noSwimmers) = swimmers.FMM.computeVel(force);

for i = 1:noSwimmers
  points(swIndex(i):swIndex(i+1)-1,:) = points(swIndex(i):swIndex(i+1)-1,:)-swimmers.x0{i}.*ones(swSizes(i),1);
end

points = reshape(points.',[],1);
force = reshape(force.',[],1);

b(1:end-6*noSwimmers) = b(1:end-6*noSwimmers) - sum(kron(ones(N,1),reshape(U,3,[])),2);

Ze  = zeros(N, noSwimmers);
x1m = zeros(N, noSwimmers);
x2m = zeros(N, noSwimmers);
x3m = zeros(N, noSwimmers);
for n = 1 : noSwimmers
    [x1,x2,x3]=extractComponents(swimmers.getIndSwimmer(n));
    x1m(swIndex(n):swIndex(n+1)-1,n) = x1;
    x2m(swIndex(n):swIndex(n+1)-1,n) = x2;
    x3m(swIndex(n):swIndex(n+1)-1,n) = x3;
end
[Omx, Omy, Omz] = extractComponents(Omega);
rotation = [ Ze -x3m x2m; x3m Ze -x1m; -x2m x1m Ze]*[Omx; Omy; Omz];
b(1:end-6*noSwimmers) = b(1:end-6*noSwimmers) + interleaveComponents(rotation(1:N),...
                                               rotation(N+1:2*N),...
                                               rotation(2*N+1:end));  

if Nearest
    af = zeros(noSwimmers, N);
    for n= 1 : noSwimmers
        af(n,swIndex(n):swIndex(n+1)-1) = ...
            sum(NN(1:Q,swIndex(n):swIndex(n+1)-1),1);
    end
    b(end-6*noSwimmers+1:end-3*noSwimmers) = kron(af,eye(3))*force;
    
    Ze  = zeros(noSwimmers,N);
    x1m = zeros(noSwimmers,N);
    x2m = zeros(noSwimmers,N);
    x3m = zeros(noSwimmers,N);
    [x1,x2,x3]=ExtractComponents(finepoints'*NNsw);
    for n=1:nSw
        x1m(n,swIndex(n):swIndex(n+1)-1)=x1(swIndex(n):swIndex(n+1)-1);
        x2m(n,swIndex(n):swIndex(n+1)-1)=x2(swIndex(n):swIndex(n+1)-1);
        x3m(n,swIndex(n):swIndex(n+1)-1)=x3(swIndex(n):swIndex(n+1)-1);
    end
    [fx, fy, fz] = extractComponents(force);
    rotation = [Ze -x3m x2m; x3m Ze -x1m; -x2m x1m Ze] * [fx;fy;fz];
    b(end-3*noSwimmers+1:end) = interleaveComponents(...
    rotation(1:noSwimmers),rotation(noSwimmers+1:2*noSwimmers),...
    rotation(2*noSwimmers+1:end));
else
    af = zeros(noSwimmers, N);
    for n= 1 : noSwimmers
        af(n,swIndex(n):swIndex(n+1)-1) = ones(swIndex(n+1)-swIndex(n),1);
    end
    b(end-6*noSwimmers+1:end-3*noSwimmers) = kron(af,eye(3))*force;

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
    [fx, fy, fz] = extractComponents(force);
    rotation = [Ze -x3m x2m; x3m Ze -x1m; -x2m x1m Ze] * [fx;fy;fz];
    b(end-3*noSwimmers+1:end) = interleaveComponents(...
    rotation(1:noSwimmers),rotation(noSwimmers+1:2*noSwimmers),...
    rotation(2*noSwimmers+1:end));
end
end

