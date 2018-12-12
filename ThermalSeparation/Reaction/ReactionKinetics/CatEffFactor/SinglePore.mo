within ThermalSeparation.Reaction.ReactionKinetics.CatEffFactor;
model SinglePore "catalyst effectiveness factor for a single pore"
extends BaseEffectiveness;
  replaceable model Thiele =
   ThieleModulus.FirstOrder constrainedby ThieleModulus.BaseThiele
                                                             annotation(choicesAllMatching=true);
    Thiele thiele(redeclare package MediumLiquid=MediumLiquid, k_v=k_v, propsLiq=propsLiq, tau=tau, epsilon=epsilon, S_x = S_x, V_cat = V_cat);
Real eta_i[nS_reac] "catalyst effectiveness factor for all reaction components";

equation
   eta = min(eta_i);
  for i in 1:nS_reac loop
    eta_i[i]=tanh(thiele.thiele[i])/thiele.thiele[i];
  end for;

end SinglePore;
