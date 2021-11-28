function [R] = Gen_Resistance(fmm)

[FU, MU] = CalculateForceMatrices(fmm);
[FO, MO] = CalculateMomentMatrices(fmm);

R = [FU, FO; MU, MO];

end

