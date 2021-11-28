function [FU, MU] = CalculateForceMatrices(fmm)

[F1, M1] = Rigid_Resistance(fmm, [1; 0; 0], [0; 0; 0]);
[F2, M2] = Rigid_Resistance(fmm, [0; 1; 0], [0; 0; 0]);
[F3, M3] = Rigid_Resistance(fmm, [0; 0; 1], [0; 0; 0]);

FU = [F1,F2,F3];
MU = [M1,M2,M3];
end

