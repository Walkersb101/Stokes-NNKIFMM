function [F, M, flag,relres,iter,resvec] = Rigid_Resistance(fmm, vel, mom, varargin)

Q = size(fmm.tree.points,1);

x = fmm.tree.points(:,1);
y = fmm.tree.points(:,2);
z = fmm.tree.points(:,3);

ux = vel(1)*ones(Q,1);
uy = vel(2)*ones(Q,1);
uz = vel(3)*ones(Q,1);

point_vel = [ux; uy; uz] + [mom(2) * z - mom(3) * y; mom(3) * x - mom(1) * z; mom(1) * y - mom(2) * x];

point_vel = interleaveComponents(point_vel(1:Q),point_vel(Q+1:2*Q),point_vel(2*Q+1:3*Q));

[f,flag,relres,iter,resvec] = fmm.computePot(point_vel,varargin{:});

if fmm.tree.arguments.nearest
    NN = kron(fmm.tree.NN,speye(3));

    [Fx, Fy, Fz] = extractComponents(NN * f);

    fine = reshape(fmm.tree.finePoints',[],1);

    x = fine' * NN;
else
    [Fx, Fy, Fz] = extractComponents(f);
    x = reshape(fmm.tree.points',[],1);
end

F = [sum(Fx);sum(Fy);sum(Fz)];
[x, y, z] = extractComponents(x);
[fx, fy, fz] = extractComponents(f);

M = [0*x, -z, y ; z, 0*x, -x; -y, x, 0*x] * [fx; fy; fz];
end