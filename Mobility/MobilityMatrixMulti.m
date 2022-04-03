function [b] = mobilityMatrixMulti(f,swimmers,B,BT)
%mobilityMatrixMulti Compute matrix vector product for the mobility problem
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
swSizes = swimmers.getSwimmerSizes();
swIndex = zeros(1,noSwimmers+1);swIndex(1) = 1;
for i = 2:noSwimmers+1
    swIndex(i) = swIndex(i-1) + swSizes(i-1);
end

force = f(1:end-6*noSwimmers);
UOm = f(end-6*noSwimmers+1:end);

N=numel(force)/3;

b = zeros(3*N+6*noSwimmers,1);

b(1:end-6*noSwimmers) = swimmers.FMM.computeVel(force);

b(1:end-6*noSwimmers) = b(1:end-6*noSwimmers) + BT*UOm;                                                                    

b(end-6*noSwimmers+1:end) = B*force;
end

