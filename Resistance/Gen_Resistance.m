function [R] = Gen_Resistance(fmm,varargin)
%Gen_Resistance Calculate the grand resistance matrix
%   Given a KIFMM object the grand resistance matrix is returned, see 
%   https://doi.org/10.1017/CBO9780511624124 for details
%
% Inputs:
%   fmm   : A KIFMM object containg all points on the body
%
% Outputs:
%   R : Grand resistance matrix for rigid body defined in creation of fmm

[F1, M1] = Rigid_Resistance(fmm, [1; 0; 0], [0; 0; 0],varargin{:});
[F2, M2] = Rigid_Resistance(fmm, [0; 1; 0], [0; 0; 0],varargin{:});
[F3, M3] = Rigid_Resistance(fmm, [0; 0; 1], [0; 0; 0],varargin{:});

FU = [F1,F2,F3];
MU = [M1,M2,M3];

[F1, M1] = Rigid_Resistance(fmm, [0; 0; 0], [1; 0; 0],varargin{:});
[F2, M2] = Rigid_Resistance(fmm, [0; 0; 0], [0; 1; 0],varargin{:});
[F3, M3] = Rigid_Resistance(fmm, [0; 0; 0], [0; 0; 1],varargin{:});

FO = [F1,F2,F3];
MO = [M1,M2,M3];

R = [FU, FO; MU, MO];
end

