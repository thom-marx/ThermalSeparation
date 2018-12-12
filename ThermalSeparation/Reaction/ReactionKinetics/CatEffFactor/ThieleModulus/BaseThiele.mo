within ThermalSeparation.Reaction.ReactionKinetics.CatEffFactor.ThieleModulus;
partial model BaseThiele "base class for thiele modulus"

final parameter Integer nS= MediumLiquid.nSubstance "number of components";
parameter SI.Volume V_cat=0.01 "Gesamtvolumen eines Katalysatorpellets";
parameter SI.Area S_x=0.02 "external catalyst surface";
replaceable package MediumLiquid = Media.BaseMediumLiquidReaction;
input Real k_v "reaction rate constant per catalyst volume";
input MediumLiquid.ThermodynamicProperties propsLiq;

    parameter Real tau "catalyst tortuosity";
    parameter Real epsilon "catalyst porosity";

output Real thiele[nS];

end BaseThiele;
