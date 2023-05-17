within ThermalSeparation.Reaction.ReactionKinetics.CatEffFactor;
model CatalystPellet "catalyst effectiveness factor for a catalyst pellet"
extends BaseEffectiveness;
Real eta_i[nS_reac] "catalyst effectiveness factor for all reaction components";

replaceable model Thiele =
    ThermalSeparation.Reaction.ReactionKinetics.CatEffFactor.ThieleModulus.FirstOrder constrainedby ThieleModulus.BaseThiele
                              annotation(choicesAllMatching=true);
    Thiele thiele(redeclare package MediumLiquid=MediumLiquid, k_v=k_v, propsLiq=propsLiq, tau=tau, epsilon=epsilon, S_x = S_x, V_cat = V_cat);

equation
eta=min(eta_i);
  for i in 1:nS_reac loop
  eta_i[i]=1/thiele.thiele[i]*(1/Modelica.Math.tanh(3*thiele.thiele[i])-1/(3*thiele.thiele[i]));
  end for;

end CatalystPellet;
