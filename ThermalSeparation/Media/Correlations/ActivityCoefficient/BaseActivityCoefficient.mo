within ThermalSeparation.Media.Correlations.ActivityCoefficient;
partial model BaseActivityCoefficient "base model for activity coefficient"
  parameter Integer nS= 2 annotation(Dialog(enable=false));
  input SI.Temperature T;
  input SI.MoleFraction x_l[nS];
  output Real gamma[nS];

end BaseActivityCoefficient;
