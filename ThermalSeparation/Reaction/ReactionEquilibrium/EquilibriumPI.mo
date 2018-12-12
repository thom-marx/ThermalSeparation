within ThermalSeparation.Reaction.ReactionEquilibrium;
model EquilibriumPI "PI controller to ensure that reaction is in equilibrium"
//homogeneous reaction
//Skript Keil, S. 126 - 148 und Bsp. S. 272

extends ThermalSeparation.Reaction.BaseReaction(film=false, redeclare
      replaceable package MediumLiquid =
        ThermalSeparation.Media.BaseMediumLiquidReaction);
 //Bsp: 2NO --> N2 + O2   nS=1--> NO, nS=2--> N2, nS=3--> O2    r_N2 = k*NO    reacComp = 2
   constant Integer nR=MediumLiquid.nR "number of reactions";
 constant Real nu[nR,nS]=MediumLiquid.nu
    "stocheometric coefficients for the reaction";
  constant Integer reacComp[nR]=MediumLiquid.reacComp
    "index of component in the component vector for which the reaction rate equation is written";
  MediumLiquid.EquilibriumConstant equilibriumConstant(c=c, gamma=gamma, propsLiq=propsLiq);
  MediumLiquid.MolarReactionEnthalpy molarReactionEnthalpy(T=T);
  Units.MolarEnthalpy h_R[nR] = molarReactionEnthalpy.h_R
    "molar reaction enthalpy, negative if exotherm";
  Real deltaH_singleReac[nR];

 Real K_eq_set[nR]= equilibriumConstant.K_eq "equilibrium constant";

  SI.MolarFlowRate NdotR[nR,nS]
    "molar flow rate for each substance and each reaction";

 /*** reaction is in equilibrium ***/
parameter Real k[nR] = fill(150000,nR) "gain for PI controller";
Real K_eq_meas[nR]= equilibriumConstant.K_eq_MWG;

   ThermalSeparation.Utilities.PI_Input[
                                nR] PI(
    each y_start=0.05,
    k=k,
   T=10)
    annotation (Placement(transformation(extent={{-28,26},{-8,46}})));

equation
  for m in 1:nR loop

    PI[m].u = (-K_eq_meas[m] + K_eq_set[m]);
    for i in 1:nS loop
         NdotR[m,i] =V*nu[m,i]/abs(nu[m,reacComp[m]])*PI[m].y;
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
end EquilibriumPI;
