function [F, M, flag,relres,iter,resvec] = Rigid_Resistance(fmm,vel,mom,...
                                                            varargin)

%Rigid_Resistance Solves the rigid body resistance problem
%   Given a KIFMM object the total force and moment are returned 
%
% Inputs:
%   fmm   : A KIFMM object containg all points on the body
%   vel   : (1,3) array containing velocity of rigid body
%   mom   : (1,3) array containing angular velocity of rigid body
%   other : Passed directly to GMRES, See gmres for details
%
% Outputs:
%   F : total force on the body
%   M : total moment on the body
%   other : See gmres for details on output                                                     
                                                        
point_vel = velocityModes(fmm.tree.points,vel,mom);

point_vel = interleaveComponents(point_vel(:,1),point_vel(:,2),...
                                 point_vel(:,3));

[f,flag,relres,iter,resvec] = gmres(@(x) this.computeVel(x),point_vel,...
                                    varargin{:});

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