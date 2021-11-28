function [FO, MO] = CalculateMomentMatrices(fmm)

[F1, M1] = Rigid_Resistance(fmm, [0; 0; 0], [1; 0; 0]);
[F2, M2] = Rigid_Resistance(fmm, [0; 0; 0], [0; 1; 0]);
[F3, M3] = Rigid_Resistance(fmm, [0; 0; 0], [0; 0; 1]);

FO = [F1,F2,F3];
MO = [M1,M2,M3];
end

