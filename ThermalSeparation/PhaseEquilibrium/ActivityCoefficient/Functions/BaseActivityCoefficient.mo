within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions;
model BaseActivityCoefficient "base model for activity coefficient"
  parameter Integer n = 1 annotation(Dialog(enable=false));
  parameter Integer nS= 2 annotation(Dialog(enable=false));
  input SI.Temperature T[n];
  input SI.MoleFraction x_l[n,nS];
  output Real gamma[n,nS];
equation

end BaseActivityCoefficient;
