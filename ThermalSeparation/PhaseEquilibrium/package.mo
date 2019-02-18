within ThermalSeparation;
package PhaseEquilibrium "Models for the Phase Equilibrium between the two Phases"
      extends Icons.Library.Orange;
  import SI = Modelica.SIunits;


  model CO2_CO2_MEA_StartUpReboiler "StartUp CO2-CO2-MEA eq. incl. reaction for Reboiler"
   extends ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium;

  Real eq_koeff[9];
  Real alpha(start=0.5); // loading in mol CO2 / mol MEA
  input Boolean startUp;
  Modelica.SIunits.Pressure p_sat_H2O(start=1.3e5);
  Modelica.SIunits.Temperature Theta = T-273.15;

     parameter Real omega_k=0.05 "large value if change between constant and variable shall be steep";
     parameter Real omega_time = 500 "Wendepunkt der tanh-Funktion";
     Real omega = 0.5*tanh(omega_k*(time-omega_time))+0.5;

  equation
  /* equilibirum of CO2 and carbamate (cf. Oexmann) */

      //if startUp then
      // x_v[mapping[2,1]]=0.2;
      //else
        Modelica.Math.log(max(1e-7,x_v[mapping[2,1]]*p)) =eq_koeff[1] + eq_koeff[2] / T + eq_koeff[3] * alpha + eq_koeff[4] * alpha / T+ eq_koeff[5] * alpha^2 +
                      eq_koeff[6] * alpha^2 / T + eq_koeff[7] * alpha^3 + eq_koeff[8] * alpha^3 / T + eq_koeff[9] * alpha^4;
      //end if;

       eq_koeff = {22.53,-7904,105.0,-16810,-286.4,26480,381.70,8295,-257.4};

    // alpha = x_l[2] / x_l[3]; // c_carba = c_CO2;

    // K[nS] = x_v[nSV]/(alpha*x_l[nSL]);

  /* equilibirum of H2O (cf. Oexmann) */

     // Modelica.Math.log(p_sat_H2O) = 73.649 - 7258.2 / T - 7.3037 * Modelica.Math.log(T) + 4.1653e-6 * T^2;
     p_sat_H2O=10^(8.07131-1730.63/((T-273.15)+233.426))*133.322;

    // K[nS-1] = p_sat_H2O/p;

  /* outputs */

    // p_bubble = 1e5; // dummy value
    p_bubble=p_sat_H2O;

    K = {factor_K[1]*p_sat_H2O/p,factor_K[2]*x_v[mapping[2,1]]/(alpha*x_l[3])};
    // K = {p_sat_H2O/p,x_v[2]/(alpha*x_l[3])};
    // K = {p_sat_H2O/p,x_v[2]/((-x_l[3]*alpha)/(2*alpha-1))};
    // K = fill(1,nS);

  end CO2_CO2_MEA_StartUpReboiler;

annotation(preferedView="info", Documentation(info="<html>
<p><h4><font color=\"#008000\">Phase Equilibrium</font></h4></p>
<p>Two phases are in equilibrium if the fugacities are the same: fiV = fiL. The vapour fugacity can be expressed in terms of the partial pressure, the fugacity coefficient and the Poynting factor. The liquid fugacity depends on the liquid composition, the saturation fugacity coefficient, the activity coefficient &gamma;i and the saturation pressure.</p>
<p>fiV = yi &middot; p &middot; &phi;iV &middot; &Phi;i </p>
<p>fiL = xi &middot; p_sat &middot; &gamma;i &middot; &phi;iL_sat </p>
<p><br/>The Poyinting factor &Phi;i is to 1, which is a good assumption if the system pressure is not too high.</p>
<p>For some substances, namely ideal gases, the liquid fugacity is rather described using Henry&apos;s law. In this case fiL = xi &middot;Hi where Hi is the Henry coefficient of the substance i. The medium model has to supply the information for which component the equilibrium is calculated using Henry&apos;s law. The medium models also supply the values for the fugacity and activity coefficients. The medium models can rely on different correlations to calculate fugacity and activity coefficients. Some of these correlations can be found in <a href=\"Modelica://ThermalSeparation.Media.Correlations\">Media.Correlations</a>. </p>
<p><br/>The VLE calculations are commonly carried out, for mixtures containing strongly polar compounds or electrolytes, with hybrid models, which use an activity coefficient model for the liquid phase and a fugacity coefficient model for the vapor phase:</p>
<p><ul>
<li>yi &middot; p &middot; &phi;iV = xi &middot; fi0 &middot; &gamma;i </li>
</ul></p>
<p>The hybrid model is the best way to represent highly non-ideal liquid mixtures at low pressures.(Lit. [1], introduction part) For systems containing dissolved gases at low pressure and at small concentrations, use Henry&apos;s law. For highly non-ideal chemical systems at high pressures, use equations of state instead of activity models(not implemented).</p>
<p><i>fi0 </i>is the reference fugacity and is commonly estimated with a saturation fugacity coefficient, the saturation pressure of the component at the system temperature and the Poynting correction for pressure (for pressures under 10 bar almost =1).</p>
<p><ul>
<li>fi0 = &phi;i_sat &middot; pi_sat &middot; &Phi;i</li>
</ul></p>
<p><br/>In case of existing supcritical substances at system temperature, the reference fugacity is replaced by Henry&apos;s law. (Lit. [2], p. 108)</p>
<p>Different classes are available to calculate the phase equilibrium: </p>
<p><a href=\"Modelica://ThermalSeparation.PhaseEquilibrium.IdealGasActivityCoeffLiquid\">IdealGasActivityCoeffLiquid</a>: fugacity coefficient of the gas phase is set to one</p>
<p><a href=\"Modelica://ThermalSeparation.PhaseEquilibrium.IdealGasIdealLiquid\">IdealGasIdealLiquid</a>: fugacity coefficients and activity coefficient are set to one </p>
<p><a href=\"Modelica://ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid\">RealGasActivityCoeffLiquid</a>: both gas and liquid phase are considered to behave non-ideal</p>
<p><a href=\"Modelica://ThermalSeparation.PhaseEquilibrium.RealGasIdealLiquid\">RealGasIdealLiquid</a>: activity coefficient is set to one </p>
<p><u><font style=\"color: #008000; \">References</font></u></p>
<p>[1] Fluid Phase Equilibria 187 (2001) 397-405</p>
<p>[2] Gmehling, Kolbe, 1992, zweite Auflage, Thermodynamik</p>
</html>"));
end PhaseEquilibrium;
