within ThermalSeparation.HeatAndMassTransfer.ThermodynamicFactorLiq;
model BaseThermodynamicFactor
  parameter Integer n = 1;
  parameter Integer nS= 2;
  constant Real R=Modelica.Constants.R;
  input SI.Temperature T[n];
  input SI.MoleFraction x[n,nS];
  output Real Gamma[n,nS,nS];
equation

end BaseThermodynamicFactor;
