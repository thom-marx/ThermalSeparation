within ThermalSeparation.Reaction;
partial model BaseReaction "partial model for reaction"

 parameter Integer n(min=1) annotation(Dialog(enable=false));
 parameter Integer nS=MediumLiquid.nSubstance
    "number of components, educts=negative, products=positive";

 input SI.Concentration c[nS];
  input SI.Volume V "liquid volume";
 input SI.MolarFlowRate Ndot_l_transfer[nS]
    "only needed if no reaction kinetic is taken into account";
  input Real gamma[nS];

  replaceable package MediumLiquid = Media.BaseMediumLiquid annotation(Dialog(enable=false));
  replaceable record Geometry = ThermalSeparation.Geometry.BasicGeometry;
 output SI.MolarFlowRate Ndot[nS]
    "molar flow rate for each substance, positive: substance is created, negative: substance is removed";
  output SI.HeatFlowRate deltaH_R "reaction enthalpy";

  // to be provided by extending class:
  constant Boolean film
    "true if reaction reaction model may also be used for film reaction";

  input MediumLiquid.ThermodynamicProperties propsLiq;
protected
  SI.MoleFraction x[nS] = propsLiq.x;
  SI.Temperature T = propsLiq.T;
  SI.Pressure p = propsLiq.p;
  annotation (Icon(graphics));
end BaseReaction;
