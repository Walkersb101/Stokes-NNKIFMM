function [F, M] = Rigid_Resistance(fmm, vel, mom)

Q = size(fmm.tree.points,1);

x = fmm.tree.points(:,1);
y = fmm.tree.points(:,2);
z = fmm.tree.points(:,3);

ux = vel(1)*ones(Q,1);
uy = vel(2)*ones(Q,1);
uz = vel(3)*ones(Q,1);

point_vel = [ux; uy; uz] + [mom(2) * z - mom(3) * y; mom(3) * x - mom(1) * z; mom(1) * y - mom(2) * x];

point_vel = Interleave_Components(point_vel(1:Q),point_vel(Q+1:2*Q),point_vel(2*Q+1:3*Q));

f = gmres(@fmm.computeVel,point_vel,10,1e-10,100);

[Fx, Fy, Fz] = Extract_Components(f);

F = [sum(Fx);sum(Fy);sum(Fz)];

M = [0*x', -z', y' ; z', 0*x', -x'; -y', x', 0*x'] * [Fx; Fy; Fz];
end