within ThermalSeparation.PhaseEquilibrium;
model H2O_CO2_MEA "H2O and CO2-MEA Eq. for CO2 separation with MEA"

 extends ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium;

Real eq_koeff[9];
Real alpha(start=0.5); // loading in mol CO2 / mol MEA
Modelica.SIunits.Pressure p_sat_H2O(start=1.3e5);
Modelica.SIunits.Temperature Theta = T-273.15;

   parameter Real omega_k=0.05
    "large value if change between constant and variable shall be steep";
   parameter Real omega_time = 500 "Wendepunkt der tanh-Funktion";
   Real omega = 0.5*tanh(omega_k*(time-omega_time))+0.5;

equation
/* equilibirum of CO2 and carbamate (cf. Oexmann) */

     Modelica.Math.log(max(1e-7,x_v[mapping[2,1]]*p)) = eq_koeff[1] + eq_koeff[2] / T + eq_koeff[3] * alpha + eq_koeff[4] * alpha / T + eq_koeff[5] * alpha^2 +
                    eq_koeff[6] * alpha^2 / T + eq_koeff[7] * alpha^3 + eq_koeff[8] * alpha^3 / T + eq_koeff[9] * alpha^4;

     eq_koeff = {22.53,-7904,105.0,-16810,-286.4,26480,381.70,8295,-257.4};

  // alpha = x_l[2] / x_l[3]; // c_carba = c_CO2;

  // K[nS] = x_v[nSV]/(alpha*x_l[nSL]);

/* equilibirum of H2O (cf. Oexmann) */

   Modelica.Math.log(p_sat_H2O) = 73.649 - 7258.2 / T - 7.3037 * Modelica.Math.log(T) + 4.1653e-6 * T^2;
   // p_sat_H2O=10^(8.07131-1730.63/((T-273.15)+233.426))*133.322;

  // K[nS-1] = p_sat_H2O/p;

/* outputs */

  p_bubble = 1e5; // dummy value

  K = {factor_K[1]*p_sat_H2O/p,factor_K[2]*x_v[mapping[2,1]]/(alpha*x_l[3])};
  // K = {p_sat_H2O/p,x_v[2]/(alpha*x_l[3])};
  // K = {p_sat_H2O/p,x_v[2]/((-x_l[3]*alpha)/(2*alpha-1))};
  // K = fill(1,nS);

end H2O_CO2_MEA;
