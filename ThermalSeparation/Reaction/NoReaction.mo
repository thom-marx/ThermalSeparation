within ThermalSeparation.Reaction;
model NoReaction "no reaction takes place"
  extends ThermalSeparation.Reaction.BaseReaction(final film=false);

equation
  Ndot = fill(0,nS);
  deltaH_R = 0;

end NoReaction;
