function [b] = mobilityPrecon(f,swimmers,B,BT,A,tol,maxIter)
%mobilityPrecon Preconditioner for Mobility matrix 
%
% Inputs:
%   f        : A (3*N + 6*noSw,1) vector of force, velocity and angular
%              velocities
%   swimmers : A Swimmer class construction containing swimmers, boundary
%              and KIFMM object
%   B        : Lower matrix needed to compute total force and torque
%   BT       : Upper matrix need to compute the velocity of individual
%              force points from the total velocity and angular velocity
%
% Outputs:
%   b : resultant vector of velocities, total force and total torque for
%       given problem

noSwimmers = swimmers.swimmerNo();

force = f(1:end-6*noSwimmers);
UOm = f(end-6*noSwimmers+1:end);

N=numel(force)/3;

b = zeros(3*N+6*noSwimmers,1);

swimmers.FMM.arguments.GMRES = 1;
[GMRESSol,~,~] = gmres(@(x) swimmers.FMM.reducedComputeVel(x), force,[],...
                        tol,maxIter);
b(1:end-6*noSwimmers) = GMRESSol;

out = -A*(B*swimmers.FMM.computeVelSparse(BT*(A*UOm)));

b(end-6*noSwimmers+1:end) = out;

[GMRESSol,~,~] = gmres(@(x) swimmers.FMM.reducedComputeVel(x), BT*out,...
                        [],tol,maxIter);
b(1:end-6*noSwimmers) = b(1:end-6*noSwimmers) - GMRESSol;
end