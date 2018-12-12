within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions;
model Wilson "multicomponent mixture after Wilson eq."
  extends BaseActivityCoefficient;

  parameter Real Lambda[nS,nS]=fill(1,nS,nS)
    "matrix with the binary coefficients";

equation
 for K in 1:n loop
  for k in 1:nS loop
      gamma[K, k] =
        ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions.WilsonFun(
        nS,
        k,
        x_l[K, :],
        Lambda);
  end for;
 end for;

end Wilson;
