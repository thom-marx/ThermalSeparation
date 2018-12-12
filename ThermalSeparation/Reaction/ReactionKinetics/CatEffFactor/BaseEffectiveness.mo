within ThermalSeparation.Reaction.ReactionKinetics.CatEffFactor;
partial model BaseEffectiveness "base class for effectiveness factor"
final parameter Integer nS_reac=MediumLiquid.nS_reac
    "number of reaction components";
replaceable package MediumLiquid = Media.BaseMediumLiquidReaction;
input Real k_v "reaction rate constant per catalyst volume";
input MediumLiquid.ThermodynamicProperties propsLiq;
/*** catalyst data ***/
parameter SI.Volume V_cat=0.01 "Gesamtvolumen eines Katalysatorpellets";
parameter SI.Area S_x=0.02 "external catalyst surface";
    parameter Real tau "catalyst tortuosity";
    parameter Real epsilon "catalyst porosity";

    output Real eta;
equation

end BaseEffectiveness;
