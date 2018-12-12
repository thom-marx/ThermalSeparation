within ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient;
model Ideal
  "fugacity coefficient for ideal liquid where the vapour phase of the single components can be regarded as ideal gas"
  import ThermalSeparation;
  extends
    ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient.BaseFugacityCoefficient;
equation
  phi_sat= ones(nS);
end Ideal;
