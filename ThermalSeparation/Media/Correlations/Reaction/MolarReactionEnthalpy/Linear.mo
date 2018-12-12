within ThermalSeparation.Media.Correlations.Reaction.MolarReactionEnthalpy;
model Linear "linear temperature dependency"
extends BaseMolarReactionEnthalpy;

parameter Real a_R[nR] =  {38.814} "h_R [J/mol] = a_R * T [K] + b_R";
parameter Real b_R[nR] = {-89155} "h_R[J/mol] = a_R * T [K] + b_R";

equation
    for m in 1:nR loop
     h_R[m] =a_R[m]*T + b_R[m];
          //if reaction is exotherm, deltaH_R0 is negative, but an energy input in the energy equation is needed
    end for;

end Linear;
