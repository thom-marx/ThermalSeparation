within ThermalSeparation.Reaction.ReactionRateCoeff;
partial model BaseReactionRateCoeff "base class for reaction rate coefficient"
  parameter Integer nR=1 annotation(Dialog(enable=false));
  input SI.Temperature T;
  output Real k[nR];
equation

end BaseReactionRateCoeff;
