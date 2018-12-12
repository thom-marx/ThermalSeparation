within ThermalSeparation.Media.Correlations.Reaction.MolarReactionEnthalpy;
model EnthalpyOfFormation
  "medium model supplies enthalpies of formation, no explicit reaction enthalpy"
extends BaseMolarReactionEnthalpy;
equation
h_R = zeros(nR);
end EnthalpyOfFormation;
