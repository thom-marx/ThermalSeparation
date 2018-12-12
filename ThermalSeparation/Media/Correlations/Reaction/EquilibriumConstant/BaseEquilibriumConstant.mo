within ThermalSeparation.Media.Correlations.Reaction.EquilibriumConstant;
partial model BaseEquilibriumConstant
  "base class for reaction equilibrium constants"
  parameter Integer nS(min=2) "number of substances";
  parameter Integer nR(min=1) "number of reactions";
  input SI.Concentration c[nS];
  input Real gamma[nS];

  input SI.MoleFraction x[nS];
  input SI.Temperature T;
  parameter SI.MolarMass MMX[nS];

  output Real K_eq[ nR] "equilibrium constant";
  output Real K_eq_MWG[ nR] "equilibrium constant";
end BaseEquilibriumConstant;
