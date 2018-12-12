within ThermalSeparation.Media.Correlations.Reaction.EquilibriumConstant;
model Activity "using activitiess for equilibrium constants"
//homogeneous reaction
//Skript Keil, S. 126 - 148 und Bsp. S. 272
//reaction kinetics taken into account

extends BaseEquilibriumConstant;
 //Bsp: 2NO --> N2 + O2   nS=1--> NO, nS=2--> N2, nS=3--> O2    r_N2 = k*NO    reacComp = 2

  parameter Real powerCoeffForward[nR,nS]
    "power coefficients for forward reaction";
 parameter Real powerCoeffReverse[nR,nS]
    "power coefficients for reverse reaction";

  Real A_forward[nR,nS];
  Real A_reverse[nR,nS];
  Real K_eq_A[ nR] = K_eq_MWG "equilibrium constant";
   Real a[nS] = x.*gamma;
equation
  for m in 1:nR loop
    for i in 1:nS loop
      A_forward[m,i] = a[i]^powerCoeffForward[m,i];
      A_reverse[m,i] = max(1e-7,a[i])^powerCoeffReverse[m,i];
    end for;
        K_eq_A[m] = product(A_reverse[m,:])/product(A_forward[m,:]);
    end for;

    //    for m in 1:nR loop
//      A[m] = -Modelica.Math.log(K0[m]) * T0[m];
//    end for;

     for m in 1:nR loop
       K_eq[m] = exp(-1/T);
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
