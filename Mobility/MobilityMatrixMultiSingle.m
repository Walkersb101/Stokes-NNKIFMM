function [b] = MobilityMatrixMulti(fmm,f,initalCon)

Nearest = fmm.tree.arguments.nearest;
points = fmm.tree.getPotentialPoints();


x0=initalCon(1:3);
b1=initalCon(4:6);
b2=initalCon(7:9);
b3=cross(b1,b2);
B=[b1(:) b2(:) b3(:)];

points = (B * points')';
points = points + x0'.*ones(size(points,1),1);

force = f(1:end-6);
U = f(end-5:end-3);
Omega = f(end-2:end);

N=numel(points)/3;

if Nearest
    finepoints = fmm.tree.getFinePoints();
    finepoints = (B * finepoints')';
    finepoints = finepoints + x0'.*ones(size(finepoints,1),1);
    finepoints = reshape(finepoints.',[],1);
end

b = zeros(3*N+6,1);

b(1:end-6) = fmm.computeVel(force);

points = reshape(points.',[],1);
force = reshape(force.',[],1);

b(1:end-6) = b(1:end-6) - kron(ones(N,1),U);

[x, y, z]=extractComponents(points);
rotation = [0*x -z y; z 0*y -x; -y x 0*z] * Omega;
b(1:end-6) = b(1:end-6) + interleaveComponents(rotation(1:N),...
                                               rotation(N+1:2*N),...
                                               rotation(2*N+1:end));

if Nearest
    NN = fmm.tree.NN;
    b(end-5:end-3)=kron(sum(NN,1),eye(3)) * force;
    [x, y, z]=extractComponents(finepoints'*kron(NN,speye(3)));
    [fx, fy, fz] = extractComponents(force);
    b(end-2:end)=[0*x -z y; z 0*y -x; -y x 0*z] * [fx;fy;fz];
else
    b(end-5:end-3) = ones(3,3*N) * force;
    [x, y, z] = extractComponents(points');
    [fx, fy, fz] = extractComponents(force);
    b(end-2:end)=[0*x -z y; z 0*y -x; -y x 0*z] * [fx;fy;fz];
end
end

