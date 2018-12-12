within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient;
model BaseFugacityCoefficient "base model for fugacity coefficients"
  parameter Integer n   annotation(Dialog(enable=false));
  parameter Integer nS annotation(Dialog(enable=false));
  input SI.Pressure p[n];
  input SI.Temperature T[n];
  input SI.MoleFraction x[n,nS];
  output Real phi[n,nS];
equation

end BaseFugacityCoefficient;
