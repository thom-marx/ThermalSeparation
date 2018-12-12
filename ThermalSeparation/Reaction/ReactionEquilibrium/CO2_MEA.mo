within ThermalSeparation.Reaction.ReactionEquilibrium;
model CO2_MEA "CO2,MEA: Heat of reaction for CO2 separation"

  extends ThermalSeparation.Reaction.BaseReaction(film=false);

// reaction is accounted for in equilibrium balance, here only heat of reaction is calculated

Real alpha; // loading in mol CO2 / mol MEA
Real h_koeff[4];
Real spec;

   parameter Real omega_k=0.05
    "large value if change between constant and variable shall be steep";
   parameter Real omega_time = 50 "Wendepunkt der tanh-Funktion";
   Real omega = 0.5*tanh(omega_k*(time-omega_time))+0.5;

constant Real R=Modelica.Constants.R;

equation
  // Ndot = fill(0,n,nS); // no bulk reaction considered - reaction only at phase boundary
//   Ndot_film = fill(0,n,nS);
//   E = fill(1,n,nS);
//   deltaH_R_film = fill(0,n);

  Ndot[:] = {0,0,0}; // MEA-Concentration is assumed constant because no carbamate is considered
  // Ndot_film[i,:] = {0,0,(-0.5*Ndot_l_transfer[i,2])}; // for each mole of carbamate created two moles of MEA are removed

/* calculation of reaction and absorption enthalpy (cf. Oexmann) */

    deltaH_R = Ndot_l_transfer[2] * (-R) * (h_koeff[1] + h_koeff[2] * alpha + h_koeff[3] *alpha^2 + h_koeff[4] * alpha^3);

    spec= (-R) * (h_koeff[1] + h_koeff[2] * alpha + h_koeff[3] *alpha^2 + h_koeff[4] * alpha^3);

    alpha = x[2]/x[3];

  h_koeff[:] = {-7904,-16810,26480,8295};

end CO2_MEA;
