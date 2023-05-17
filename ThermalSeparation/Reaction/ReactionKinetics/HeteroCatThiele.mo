within ThermalSeparation.Reaction.ReactionKinetics;
model HeteroCatThiele "heterogeneous catalysis including thiele modulus"

extends BaseReaction(redeclare replaceable package MediumLiquid =
        ThermalSeparation.Media.BaseMediumLiquidReaction, film=false);
  constant Integer nR=MediumLiquid.nR "number of reactions";
  constant Real nu[nR,nS]=MediumLiquid.nu
    "stocheometric coefficients for the reaction";
  constant Integer reacComp[nR]=MediumLiquid.reacComp
    "index of component in the component vector for which the reaction rate equation is written";
  MediumLiquid.ReactionRates reactionRates(c=c, gamma=gamma, propsLiq=propsLiq);
  SI.MolarFlowRate NdotR[nR,nS]
    "molar flow rate for each substance and each reaction";
  ThermalSeparation.Units.ReactionRate r_intrinsic[nR]=m_cat/n*reactionRates.r
    "intrinsic reaction rate";
  ThermalSeparation.Units.ReactionRate r[nR];
  Real eta=catEff.eta "catalyst effectiveness factor";

  /*** reaction enthalpy ***/
   MediumLiquid.MolarReactionEnthalpy molarReactionEnthalpy(T=T);
   Units.MolarEnthalpy h_R[nR] = molarReactionEnthalpy.h_R
    "molar reaction enthalpy, negative if exotherm";
  Real deltaH_singleReac[nR];

  /*** catalyst data ***/
  parameter SI.Volume V_cat=0.01 "Gesamtvolumen eines Katalysatorpellets";
parameter SI.Area S_x=0.02 "external catalyst surface";
    parameter Real tau "catalyst tortuosity";
    parameter Real epsilon "catalyst porosity";
  parameter SI.Mass m_cat "mass of catalysator in the whole column in Gramm";

  replaceable model CatalystEffectiveness =
      ThermalSeparation.Reaction.ReactionKinetics.CatEffFactor.BaseEffectiveness                        annotation(choicesAllMatching=true);
  CatalystEffectiveness catEff(redeclare package MediumLiquid=MediumLiquid, k_v=reactionRates.k_forward[1], propsLiq=propsLiq, tau=tau, epsilon=epsilon, S_x = S_x, V_cat = V_cat);

equation
  for m in 1:nR loop
r[m] = r_intrinsic[m] * eta;
   for i in 1:nS loop
 NdotR[m,i] =nu[m,i]/abs(nu[m,reacComp[m]])*r[m];
  end for;
end for;
  for i in 1:nS loop
    Ndot[i] = sum(NdotR[:,i]);
  end for;

/*** reaction enthalpy ***/

    for m in 1:nR loop
    deltaH_singleReac[m] =- abs(NdotR[m,reacComp[m]])*h_R[m]/abs(nu[m,reacComp[m]]);
    end for;
 deltaH_R = sum(deltaH_singleReac[:]);
  annotation (Documentation(info="<html>
<p>This model is intended for heterogeneous catalytic reactions where the change in reaction kinetics due to diffusion in the catalysator is taken into account using the Thiele modulus. The (intrinsic) reaction rate is determined at liquid bulk conditions using the equations provided by the medium model. Using the Thiele modulus a catalyst effectiveness factor is determined which is used to determine the (effective) reaction rate.</p>
</html>"));
end HeteroCatThiele;
