within ThermalSeparation.PhaseEquilibrium.HenrysLaw;
model BaseHenry
  replaceable package MediumVapour = 
      ThermalSeparation.Media.IdealGasMixtures.N2_H2O;
  parameter Integer n;
  parameter Integer nSV;
  parameter Integer nSL;

  parameter Boolean henry_temp;

  input SI.MoleFraction x_l[n,nSL];
  input SI.Temperature T[n];
  output Real He[n,nSV];

end BaseHenry;
