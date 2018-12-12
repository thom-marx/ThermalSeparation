within ThermalSeparation.Reaction.ReactionRateCoeff;
model Arrhenius
  extends BaseReactionRateCoeff;
  parameter Real k_0[nR] = 1*ones(nR);
  parameter SI.Energy E[nR] = 2000*ones(nR);
equation

  for r in 1:nR loop
  k[r] = k_0*exp(-E/(Modelica.Constants.R*T));
  end for;

end Arrhenius;
