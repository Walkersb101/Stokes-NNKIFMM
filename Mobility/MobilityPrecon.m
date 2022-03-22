function [b] = MobilityPrecon(f,swimmers)
points = swimmers.getSwimmers();

noSwimmers = swimmers.swimmerNo();
swSizes = swimmers.getSwimmerSizes();
swIndex = zeros(1,noSwimmers+1);swIndex(1) = 1;
for i = 2:noSwimmers+1
    swIndex(i) = swIndex(i-1) + swSizes(i-1);
end

force = f(1:end-6*noSwimmers);
U = f(end-6*noSwimmers+1:end-3*noSwimmers);
Omega = f(end-3*noSwimmers+1:end);

points = reshape(points.',[],1);
force = reshape(force.',[],1);

N=numel(points)/3;

b = zeros(3*N+6*noSwimmers,1);

swimmers.FMM.arguments.GMRES = 1;
b(1:end-6*noSwimmers) = gmres(@(x) swimmers.FMM.reducedComputeVel(x), force,[],1e-4,100);

af = zeros(noSwimmers, N);
for n= 1 : noSwimmers
    af(n,swIndex(n):swIndex(n+1)-1) = -ones(swIndex(n+1)-swIndex(n),1);
end

Ze  = zeros(noSwimmers,N);
x1m = zeros(noSwimmers,N);
x2m = zeros(noSwimmers,N);
x3m = zeros(noSwimmers,N);
for n=1:noSwimmers
    [x1,x2,x3]=extractComponents(swimmers.getIndSwimmer(n));
    x1m(n,swIndex(n):swIndex(n+1)-1)=x1;
    x2m(n,swIndex(n):swIndex(n+1)-1)=x2;
    x3m(n,swIndex(n):swIndex(n+1)-1)=x3;
end
B = vertcat(kron(eye(3),af),[Ze x3m -x2m; -x3m Ze x1m; x2m -x1m Ze]);
A = B*B.';
[Ux,Uy,Uz] = extractComponents(U);
[Omx, Omy, Omz] = extractComponents(Omega);

RHP = B.'*(A\[Ux;Uy;Uz;Omx;Omy;Omz]);
RHP = interleaveComponents(RHP(1:N),RHP(N+1:2*N),RHP(2*N+1:end));

[fx,fy,fz] = extractComponents(swimmers.FMM.computeVel(RHP));
out = -A\(B*[fx;fy;fz]);

topSlice = out(end-6*noSwimmers+1:end-3*noSwimmers);
bottomSlice = out(end-3*noSwimmers+1:end);


b(end-6*noSwimmers+1:end)=[interleaveComponents(...
                           topSlice(1:noSwimmers),...
                           topSlice(noSwimmers+1:2*noSwimmers),...
                           topSlice(2*noSwimmers+1:end));...
                           interleaveComponents(...
                           bottomSlice(1:noSwimmers),...
                           bottomSlice(noSwimmers+1:2*noSwimmers),...
                           bottomSlice(2*noSwimmers+1:end))];


RHP = B.'*out;
RHP = interleaveComponents(RHP(1:N),RHP(N+1:2*N),RHP(2*N+1:end));
b(1:end-6*noSwimmers) = b(1:end-6*noSwimmers) - ...
                        gmres(@(x) swimmers.FMM.reducedComputeVel(x),...
                              RHP ,[],1e-4,100);
end

