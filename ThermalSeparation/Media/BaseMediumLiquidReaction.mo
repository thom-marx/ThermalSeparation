within ThermalSeparation.Media;
package BaseMediumLiquidReaction
  "base class for liquid medium with reacting substances"
  extends BaseMediumLiquid;
   constant Integer nR(min=1)=1 "number of reactions";

  //Bsp: A + 2B --> 0.5 C   ---> nu = {-1,-2,0.5}
 constant Real nu[nR,nSubstance]=fill(1,nR,nSubstance) "stocheometric coefficients for the reaction";

  constant Integer reacComp[nR](min=fill(1,nR))=fill(1,nR)  "index of component in the component vector for which the reaction rate equation is written, to which the molar reaction enthalpy 
    refers ... ";

constant Integer nS_reac(min=1)=1 "number of reaction components";

replaceable partial model MolarReactionEnthalpy
    "base class for molar reaction enthalpy"

  input SI.Temperature T;

  output Units.MolarEnthalpy h_R[nR]
      "molar reaction enthalpy, negative if exotherm";

end MolarReactionEnthalpy;

  replaceable partial model ReactionRates "models for reaction rates"
    input SI.Concentration c[nSubstance];
    input Real gamma[nSubstance];
    input ThermodynamicProperties propsLiq;
    output ThermalSeparation.Units.ReactionRate r[nR];
    output Real k_forward[nR] "reaction rate coefficient forward reaction";

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end ReactionRates;

  replaceable partial model EquilibriumConstant
    input SI.Concentration c[nSubstance];
    input Real gamma[nSubstance];
    input ThermodynamicProperties propsLiq;
    output Real K_eq[ nR] "reaction equilibrium constant";
    output Real K_eq_MWG[ nR]
      "reaction equilibrium constant on the basis of concentrations or activities";
  equation

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end EquilibriumConstant;
end BaseMediumLiquidReaction;
