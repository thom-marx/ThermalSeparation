within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient;
model Ideal "composition equals activity (gamma = 1)"
  extends BaseActivityCoefficient;

equation
gamma = fill(1,nS);//ones(n,nS);

end Ideal;
