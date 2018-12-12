within ThermalSeparation.Media.Correlations.Reaction;
package ReactionRates "models for reaction rates"
  partial model BaseReactionRates "base class for reaction rates"
    parameter Integer nS(min=2) "number of substances";
    parameter Integer nR(min=1) "number of reactions";
    input SI.Concentration c[nS];
    input Real gamma[nS];
   output ThermalSeparation.Units.ReactionRate r[nR];
    output Real k_forward[nR] "reaction rate coefficient forward reaction";

    input SI.MoleFraction x[nS];
    input SI.Temperature T;
    parameter SI.MolarMass MMX[nS];

  end BaseReactionRates;

model Concentration "using concentrations for reaction rate"
//homogeneous reaction
//Skript Keil, S. 126 - 148 und Bsp. S. 272
//reaction kinetics taken into account

extends BaseReactionRates;
 //Bsp: 2NO --> N2 + O2   nS=1--> NO, nS=2--> N2, nS=3--> O2    r_N2 = k*NO    reacComp = 2

  parameter Boolean reverse[nR]= fill(false,nR)
      "true if reverse reaction is taken into account";
  parameter Real powerCoeffForward[nR,nS]
      "power coefficients for forward reaction";
 parameter Real powerCoeffReverse[nR,nS]
      "power coefficients for reverse reaction";

  Real C_forward[nR,nS];
  Real C_reverse[nR,nS];
  ThermalSeparation.Units.ReactionRate r_forward[
                                    nR] "reaction rate for forward reaction";
  ThermalSeparation.Units.ReactionRate r_reverse[
                                    nR] "reaction rate for reverse reaction";

   /*** reaction rate coefficient ***/
  replaceable model RRC_forward =
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.Arrhenius
                                                             constrainedby
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.BaseReactionRateCoeff
                                                                                                        annotation(choicesAllMatching=true);
  RRC_forward rrc_forward(nR=nR, T=T);

 replaceable model RRC_reverse =
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.Arrhenius
                                                             constrainedby
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.BaseReactionRateCoeff
                                                                                                        annotation(choicesAllMatching=true);
  RRC_reverse rrc_reverse(nR=nR, T = T);
  Real k_reverse[nR]=rrc_reverse.k "reaction rate coefficient reverse reaction";

equation
     k_forward= rrc_forward.k "reaction rate coefficient forward reaction";
  for m in 1:nR loop
    for i in 1:nS loop
      C_forward[m,i] = c[i]^powerCoeffForward[m,i];
      C_reverse[m,i] = max(1e-7,c[i])^powerCoeffReverse[m,i];
    end for;
    r_forward[m] = k_forward[m]*product(C_forward[m,:]);
    r_reverse[m] = r_forward[m] - k_reverse[m]*product(C_reverse[m,:]);
   r[m] = if reverse[m] then r_reverse[m] else r_forward[m];
   end for;

  annotation (Documentation(info="<html>
<p>This model takes into account reaction. It is possible to have several (independent, parallel or consecutive) reactions, irreversible as well as reversible reactions. A kinetic approach is used, i.e. the equilibrium is not attained at once.</p>
<p>The parameter <i><b>nR</b></i> denotes the number of reactions. For each reaction, it has to be stated whether the reaction is reversible or irreversible using the boolean parameter vector <i><b>reversible</b></i>. The stocheometric coefficients of the reaction are stored in the parameter <i><b>nu</b></i>. The reaction order for forward and reverse reaction is stored in<i><b> powerCoeffForward</b></i> and <i><b>powerCoeffReverse</b></i> respectively. If no reverse reaction takes place, any values for <i>powerCoeffReverse</i> can be used. The reaction order can not vary during simulation (even though in reality there are reactions where the reaction order changes; often due to a change in temperature). The reaction rate equation is usually written for one component. The number of this component has to be stated in the parameter <i><b>reacComp</b></i>.</p>
<p>The reaction rate coefficient of the reverse reaction, <i>k_reverse</i>, is calculated using the reaction rate coefficient of the forward reaction, <i>k_forward</i>, and the equilibrium constant <i>K</i>. The values for <i>k_forward</i> and <i>K</i> are calculated in separate classes which are declared replaceable in this model. The assumption that K is the ratio of k_forward to k_reverse only holds for elementary reactions, not for overall reactions (so it is true for example for A+B &LT;--&GT; C and C+B &LT;--&GT; D but not necessarily for A+2B &LT;--&GT; D).</p>
<p><u><b>Example:</b></u></p>
<p><ul>
<li>Medium with 5 substances: A, B, C, D, E</li>
<li>reaction 1: </li>
<p>A + D &LT;--&GT; C </p>
<p>r_A = k1_forward * c_D*c_A - k1_reverse * c_C</p>
<li>reaction 2: </li>
</ul></p>
<p>B--&GT; E </p>
<p>r_B = k2 * c_B^2</p>
<p><h4>Parameters:</h4></p>
<p><ul>
<li>nR = 2</li>
<li>reversible = {true, false}</li>
<li>nu = {{-1, 0, 1, -1, 0},{0, -1, 0, 0, 1}}</li>
<li>powerCoeffForward = {{1, 0, 0, 1, 0},{0, 2, 0, 0, 0}}</li>
<li>powerCoeffReverse = {{0, 0, 1, 0, 0},{0, 0, 0, 0, 0}}</li>
<li>reacComp = {1, 2}</li>
</ul></p>
</html>"));
end Concentration;

model Activity "using activitiess for reaction rate"
//homogeneous reaction
//Skript Keil, S. 126 - 148 und Bsp. S. 272
//reaction kinetics taken into account

extends BaseReactionRates;
 //Bsp: 2NO --> N2 + O2   nS=1--> NO, nS=2--> N2, nS=3--> O2    r_N2 = k*NO    reacComp = 2

  parameter Boolean reverse[nR]= fill(false,nR)
      "true if reverse reaction is taken into account";
  parameter Real powerCoeffForward[nR,nS]
      "power coefficients for forward reaction";
 parameter Real powerCoeffReverse[nR,nS]
      "power coefficients for reverse reaction";

  Real A_forward[nR,nS];
  Real A_reverse[nR,nS];
  ThermalSeparation.Units.ReactionRate r_forward[
                                    nR] "reaction rate for forward reaction";
  ThermalSeparation.Units.ReactionRate r_reverse[
                                    nR] "reaction rate for reverse reaction";

   /*** reaction rate coefficient ***/
  replaceable model RRC_forward =
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.Arrhenius
                                                             constrainedby
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.BaseReactionRateCoeff
                                                                                                        annotation(choicesAllMatching=true);
  RRC_forward rrc_forward(nR=nR, T=T);

 replaceable model RRC_reverse =
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.Arrhenius
                                                             constrainedby
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.BaseReactionRateCoeff
                                                                                                        annotation(choicesAllMatching=true);
  RRC_reverse rrc_reverse(nR=nR, T = T);

  Real k_reverse[nR]=rrc_reverse.k "reaction rate coefficient reverse reaction";
   Real a[nS] = x.*gamma;
equation
    k_forward= rrc_forward.k "reaction rate coefficient forward reaction";
  for m in 1:nR loop
    for i in 1:nS loop
      A_forward[m,i] = a[i]^powerCoeffForward[m,i];
      A_reverse[m,i] = max(1e-7,a[i])^powerCoeffReverse[m,i];
    end for;
    r_forward[m] = k_forward[m]*product(A_forward[m,:]);
    r_reverse[m] = r_forward[m] - k_reverse[m]*product(A_reverse[m,:]);
   r[m] = if reverse[m] then r_reverse[m] else r_forward[m];
   end for;

  annotation (Documentation(info="<html>
<p>This model takes into account reaction. It is possible to have several (independent, parallel or consecutive) reactions, irreversible as well as reversible reactions. A kinetic approach is used, i.e. the equilibrium is not attained at once.</p>
<p>The parameter <i><b>nR</b></i> denotes the number of reactions. For each reaction, it has to be stated whether the reaction is reversible or irreversible using the boolean parameter vector <i><b>reversible</b></i>. The stocheometric coefficients of the reaction are stored in the parameter <i><b>nu</b></i>. The reaction order for forward and reverse reaction is stored in<i><b> powerCoeffForward</b></i> and <i><b>powerCoeffReverse</b></i> respectively. If no reverse reaction takes place, any values for <i>powerCoeffReverse</i> can be used. The reaction order can not vary during simulation (even though in reality there are reactions where the reaction order changes; often due to a change in temperature). The reaction rate equation is usually written for one component. The number of this component has to be stated in the parameter <i><b>reacComp</b></i>.</p>
<p>The reaction rate coefficient of the reverse reaction, <i>k_reverse</i>, is calculated using the reaction rate coefficient of the forward reaction, <i>k_forward</i>, and the equilibrium constant <i>K</i>. The values for <i>k_forward</i> and <i>K</i> are calculated in separate classes which are declared replaceable in this model. The assumption that K is the ratio of k_forward to k_reverse only holds for elementary reactions, not for overall reactions (so it is true for example for A+B &LT;--&GT; C and C+B &LT;--&GT; D but not necessarily for A+2B &LT;--&GT; D).</p>
<p><u><b>Example:</b></u></p>
<p><ul>
<li>Medium with 5 substances: A, B, C, D, E</li>
<li>reaction 1: </li>
<p>A + D &LT;--&GT; C </p>
<p>r_A = k1_forward * c_D*c_A - k1_reverse * c_C</p>
<li>reaction 2: </li>
</ul></p>
<p>B--&GT; E </p>
<p>r_B = k2 * c_B^2</p>
<p><h4>Parameters:</h4></p>
<p><ul>
<li>nR = 2</li>
<li>reversible = {true, false}</li>
<li>nu = {{-1, 0, 1, -1, 0},{0, -1, 0, 0, 1}}</li>
<li>powerCoeffForward = {{1, 0, 0, 1, 0},{0, 2, 0, 0, 0}}</li>
<li>powerCoeffReverse = {{0, 0, 1, 0, 0},{0, 0, 0, 0, 0}}</li>
<li>reacComp = {1, 2}</li>
</ul></p>
</html>"));
end Activity;

  model Langmuir
    "reaction rate using Langmuir sorption isotherms and activities"
  //homogeneous reaction
  //Skript Keil, S. 126 - 148 und Bsp. S. 272
  //reaction kinetics taken into account
  extends BaseReactionRates;
  //r is the reaction rate per gramm of catalyst

    parameter Real powerCoeffForward[nR,nS]
      "power coefficients for forward reaction";
   parameter Real powerCoeffReverse[nR,nS]
      "power coefficients for reverse reaction";
      parameter Real K_ads[nS];

    Real A_forward[nR,nS];
    Real A_reverse[nR,nS];
    ThermalSeparation.Units.ReactionRate r_forward[
                                      nR] "reaction rate for forward reaction";
    ThermalSeparation.Units.ReactionRate r_reverse[
                                      nR] "reaction rate for reverse reaction";
   ThermalSeparation.Units.ReactionRate r_back[
                                      nR] "reaction rate for reverse reaction";

    Real a_dash[nR,nS];
    final parameter SI.MolarMass MM[nS] = MMX * 1000;

     Real sum_a_dash[nR];

        /*** reaction rate coefficient ***/
    replaceable model RRC_forward =
        ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.Arrhenius
                                                               constrainedby
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.BaseReactionRateCoeff
                                                                                                          annotation(choicesAllMatching=true);
    RRC_forward rrc_forward(nR=nR, T=T);

   replaceable model RRC_reverse =
        ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.Arrhenius
                                                               constrainedby
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.BaseReactionRateCoeff
                                                                                                          annotation(choicesAllMatching=true);
    RRC_reverse rrc_reverse(nR=nR, T = T);

    Real k_reverse[nR]=rrc_reverse.k
      "reaction rate coefficient reverse reaction";

  equation
      k_forward= rrc_forward.k "reaction rate coefficient forward reaction";
      for m in 1:nR loop
      for i in 1:nS loop
        a_dash[m,i] = K_ads[i]/MM[i]*x[i]*gamma[i];
        A_forward[m,i] = max(1e-9,a_dash[m,i])^powerCoeffForward[m,i];
        A_reverse[m,i] = max(1e-9,a_dash[m,i])^powerCoeffReverse[m,i];
      end for;
      r_forward[m] = k_forward[m]*product(A_forward[m,:]);
     r_reverse[m] = r_forward[m] - k_reverse[m]*product(A_reverse[m,:]);
      r_back[m] = k_reverse[m]*product(A_reverse[m,:]);
      sum_a_dash[m]=sum(a_dash[m,:]);
      r[m] =   r_reverse[m]/(sum_a_dash[m]^2);

    end for;

    annotation (Documentation(info="<html>
<p>This model takes into account reaction. It is possible to have several (independent, parallel or consecutive) reactions, irreversible as well as reversible reactions. A kinetic approach is used, i.e. the equilibrium is not attained at once.</p>
<p>The parameter <i><b>nR</b></i> denotes the number of reactions. For each reaction, it has to be stated whether the reaction is reversible or irreversible using the boolean parameter vector <i><b>reversible</b></i>. The stocheometric coefficients of the reaction are stored in the parameter <i><b>nu</b></i>. The reaction order for forward and reverse reaction is stored in<i><b> powerCoeffForward</b></i> and <i><b>powerCoeffReverse</b></i> respectively. If no reverse reaction takes place, any values for <i>powerCoeffReverse</i> can be used. The reaction order can not vary during simulation (even though in reality there are reactions where the reaction order changes; often due to a change in temperature). The reaction rate equation is usually written for one component. The number of this component has to be stated in the parameter <i><b>reacComp</b></i>.</p>
<p>The reaction rate coefficient of the reverse reaction, <i>k_reverse</i>, is calculated using the reaction rate coefficient of the forward reaction, <i>k_forward</i>, and the equilibrium constant <i>K</i>. The values for <i>k_forward</i> and <i>K</i> are calculated in separate classes which are declared replaceable in this model. The assumption that K is the ratio of k_forward to k_reverse only holds for elementary reactions, not for overall reactions (so it is true for example for A+B &LT;--&GT; C and C+B &LT;--&GT; D but not necessarily for A+2B &LT;--&GT; D).</p>
<p><u><b>Example:</b></u></p>
<p><ul>
<li>Medium with 5 substances: A, B, C, D, E</li>
<li>reaction 1: </li>
<p>A + D &LT;--&GT; C </p>
<p>r_A = k1_forward * c_D*c_A - k1_reverse * c_C</p>
<li>reaction 2: </li>
</ul></p>
<p>B--&GT; E </p>
<p>r_B = k2 * c_B^2</p>
<p><h4>Parameters:</h4></p>
<p><ul>
<li>nR = 2</li>
<li>reversible = {true, false}</li>
<li>nu = {{-1, 0, 1, -1, 0},{0, -1, 0, 0, 1}}</li>
<li>powerCoeffForward = {{1, 0, 0, 1, 0},{0, 2, 0, 0, 0}}</li>
<li>powerCoeffReverse = {{0, 0, 1, 0, 0},{0, 0, 0, 0, 0}}</li>
<li>reacComp = {1, 2}</li>
</ul></p>
</html>"));
  end Langmuir;

  package ReactionRateCoeff

    model Arrhenius
      extends BaseReactionRateCoeff;
    //    parameter Real k_0[nR] = 1*ones(nR);
    //    parameter SI.Energy E[nR] = 2000*ones(nR);
    equation

      for r in 1:nR loop
      k[r] = k_0*exp(-E/(Modelica.Constants.R*T));
      end for;

    end Arrhenius;

    partial model BaseReactionRateCoeff
      "base class for reaction rate coefficient"
      parameter Integer nR=1 annotation(Dialog(enable=false));

         parameter Real k_0[nR] = 1*ones(nR);
       parameter SI.Energy E[nR] = 2000*ones(nR);
      input SI.Temperature T;
      output Real k[nR];
    equation

    end BaseReactionRateCoeff;
  end ReactionRateCoeff;
end ReactionRates;
