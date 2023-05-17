within ThermalSeparation.Reaction.ReactionKinetics;
model HeteroCatMasstransfer
  "heterogeneous catalysis including liquid-solid mass transfer"
//homogeneous reaction
//Skript Keil, S. 126 - 148 und Bsp. S. 272
//reaction kinetics taken into account
extends BaseReaction(redeclare replaceable package MediumLiquid =
        ThermalSeparation.Media.BaseMediumLiquidReaction, film=false, redeclare replaceable record Geometry =
        ThermalSeparation.Geometry.BasicGeoPackedColumn);
  constant Integer nR=MediumLiquid.nR "number of reactions";
  constant Real nu[nR,nS]=MediumLiquid.nu
    "stocheometric coefficients for the reaction";
  constant Integer reacComp[nR]=MediumLiquid.reacComp
    "index of component in the component vector for which the reaction rate equation is written";
  MediumLiquid.ReactionRates reactionRates(c=mediumLiquidCat.c, gamma=activityCoeffCat.gamma, propsLiq=mediumLiquidCat.properties);
  SI.MolarFlowRate NdotR[nR,nS]
    "molar flow rate for each substance and each reaction";
  ThermalSeparation.Units.ReactionRate r[nR]=m_cat/n*reactionRates.r;

  /*** reaction enthalpy ***/
   MediumLiquid.MolarReactionEnthalpy molarReactionEnthalpy(T=mediumLiquidCat.T);
   Units.MolarEnthalpy h_R[nR] = molarReactionEnthalpy.h_R
    "molar reaction enthalpy, negative if exotherm";
  Real deltaH_singleReac[nR];
  parameter SI.Mass m_cat "mass of catalysator in the whole column in Gramm";

/*** values on catalyst surface ***/
  SI.MoleFraction x_cat[nS] "mole fraction on catalyst surface";
  MediumLiquid.ActivityCoefficient activityCoeffCat(T=mediumLiquidCat.T,x_l=x_cat);
  ThermalSeparation.FilmModel.BaseClasses.MaxwellStefanMatrix
    maxwellStefanMatrixLiq(                                                               nS = nS, x=x, dummy=k_l);
  Real R[nS-1,nS-1] = maxwellStefanMatrixLiq.matrix;
  final parameter Integer nSi=if sum(abs(MediumLiquid.ic))==0 then nS-1 else nS-2
    "number of independent composition gradients";
  Real pot_diff;// difference in electrostatic potential
  MediumLiquid.ThermodynamicFactor thermoFactor( T=T, x=x);
  MediumLiquid.BaseProperties mediumLiquidCat(T=T, x=x_cat, p=p);
  SI.Area A=geometry.a*geometry.H*geometry.A/n "interfacial liquid-solid area";
  Geometry geometry;
  parameter Real k_l[nS,nS] = fill(1e-4,nS,nS);

equation
  for m in 1:nR loop

   for i in 1:nS loop
 NdotR[m,i] =nu[m,i]/abs(nu[m,reacComp[m]])*r[m];
  end for;
end for;
  for i in 1:nS loop
    Ndot[i] = sum(NdotR[:,i]);
  end for;

  /*** mass transfer ***/
  R[1:nSi,1:nSi]*(Ndot[1:nSi]) =sum(Ndot[:]) * R[1:nSi,1:nSi]* x[1:nSi] + A * (1/propsLiq.v* thermoFactor.Gamma[1:nSi,1:nSi]*(x_cat[1:nSi] - x[1:nSi]) + Modelica.Constants.F/(Modelica.Constants.R*T)*propsLiq.c[1:nSi].*MediumLiquid.ic[1:nSi]*pot_diff);
  sum(x_cat) = 1;
  if sum(abs(MediumLiquid.ic))==0 then
        //no charged particles are in the solution
        pot_diff =0;
  else  sum(MediumLiquid.ic[:].*x_cat[:])=0; //electroneutrality condition
    pot_diff= -sum(MediumLiquid.ic[:].*(-mediumLiquidCat.c + propsLiq.c))*(Modelica.Constants.R*T)/sum(MediumLiquid.ic[:].*MediumLiquid.ic[:].*propsLiq.c)/Modelica.Constants.F;
      end if;

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
end HeteroCatMasstransfer;
