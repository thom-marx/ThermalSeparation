within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions;
model UNIQUAC "multicomponent mixture after UNIQUAC eq."
  extends BaseActivityCoefficient;

  parameter Real q[nS]=fill(1,nS);
  parameter Real r[nS]=fill(3,nS);
  parameter Real u[nS,nS]=fill(2,nS,nS);

equation
 for K in 1:n loop
  for k in 1:nS loop
      gamma[K, k] =
        ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions.UNIQUACFun(
        nS,
        k,
        x_l[K, :],
        T[K],
        r,
        q,
        u);
  end for;
end for;

end UNIQUAC;
