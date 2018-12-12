within ThermalSeparation.Reaction.Homogeneous;
model ReactionBulkActivity "Methylacetate synthesis"
//homogeneous reaction
//Skript Keil, S. 126 - 148 und Bsp. S. 272
//reaction kinetics taken into account

extends ThermalSeparation.Reaction.BaseReaction;

  parameter Integer nR=2 "number of reactions";
//Bsp: A + 2B --> 0.5 C   ---> nu = {-1,-2,0.5}
  parameter Real nu[nR,nS]={{-2,-2,0.5},{-1,-2,0.5}}
    "stocheometric coefficients for the reaction";
  parameter Boolean reverse[nR] = fill(false,nR)
    "true if reverse reaction is taken into account";

  parameter Real powerCoeffForward[nR,nS]={{1,1,0}, {1,1,0}}
    "power coefficients for forward reaction";
  parameter Real powerCoeffReverse[nR,nS]
    "power coefficients for reverse reaction";
  parameter Integer reacComp[nR]={1,2}
    "index of component in the component vector for which the reaction rate equation is written";
//Bsp: 2NO --> N2 + O2   nS=1--> NO, nS=2--> N2, nS=3--> O2    r_N2 = k*NO    reacComp = 2
parameter SI.Energy deltaH_R0[nR] "negative if reaction is exotherm";
  Units.MolarEnthalpy h_R[nR] "molar reaction enthalpy, negative if exotherm";
  parameter Real a_R[nR] = {38.814} "h_R [J/mol] = a_R * T [K] + b_R";
  parameter Real b_R[nR] = {-89155} "h_R [J/mol] = a_R * T [K] + b_R";

//vielleicht kann man hier was mit conditional model machen (so dass man für jede Reaktion ein anderes ReactionRateCoeff-Model verwenden kann
//und überhaupt sollte man sich nochmal angucken, ob man die Temperatur so übergeben kann, oder ob man eine for-Schleife braucht
  replaceable model RRC = 
      ThermalSeparation.Reaction.ReactionRateCoeff.Arrhenius;
  RRC rrc( nR=nR, T = T, k_0=8.497e6*ones(nR), E=60.47e3*ones(nR));
  RRC rrc_back( nR=nR, T = T, k_0=6.127e5*ones(nR), E=63.73e3*ones(nR));

  replaceable model EqConst = 
      ThermalSeparation.Reaction.EquilibriumConstant.Constant constrainedby
    ThermalSeparation.Reaction.EquilibriumConstant.BaseK 
                              annotation(choicesAllMatching=true);
  EqConst eqConst( nR=nR, T=T);

protected
  Real C_forward[nR,nS];
  Real C_reverse[nR,nS];
  ThermalSeparation.Units.ReactionRate r_forward[
                                    nR] "reaction rate for forward reaction";
  ThermalSeparation.Units.ReactionRate r_reverse[
                                    nR] "reaction rate for reverse reaction";
 ThermalSeparation.Units.ReactionRate r_back[
                                    nR] "reaction rate for reverse reaction";
  Real k_forward[nR]= rrc.k "reaction rate coefficient forward reaction";
  Real k_reverse[nR]= rrc_back.k "reaction rate coefficient reverse reaction";
  Real K_eq[nR]= eqConst.K "equilibrium constant";
  ThermalSeparation.Units.ReactionRate r[
                            nR];
  SI.MolarFlowRate NdotR[        nR,nS]
    "molar flow rate for each substance and each reaction";
Real test_pCf[nR];
Real test_pCr[nR];
SI.Energy deltaH_singleReac[nR];

 /*** reaction is in equilibrium ***/
public
parameter Boolean equilibrium[nR] = fill(false,nR)
    "reaction is in equilibrium, no kinetics";
parameter Real k[nR] = fill(500,nR) "gain for PI controller";
Real K_bulk[nR];

  Modelica.Blocks.Continuous.PI[nR] PI(
    each initType=Modelica.Blocks.Types.Init.InitialOutput,
    each y_start=0.05,
    each T=10,
    k=k) 
    annotation (Placement(transformation(extent={{-28,26},{-8,46}})));

parameter SI.Mass m_cat "mass of catalysator in the whole column in Gramm";
Real a_dash[nR,nS];
//parameter SI.MolarMass MM[nS]={0.0741,0.0601,0.03204,0.018}*1000;
final parameter SI.MolarMass MM[nS] = MediumLiquid.MMX * 1000;
parameter Real K_ads[nS]={4.15, 3.15, 5.64, 5.24}; //Forner
//parameter Real K_ads[nS]={0.125, 0.24, 1.57, 1.6};

     parameter Real omega_k=0.05;// Modelica-Medien: 0.05
   // "large value if change between constant and variable shall be steep";
   parameter Real omega_time = 500;// Modelica-Medien: 500 "Wendepunkt der tanh-Funktion";
   Real omega= 1;// 0.5*tanh(omega_k*(time-omega_time))+0.5;

   Real sum_a_dash[nR];
equation

  for m in 1:nR loop
    for i in 1:nS loop
      a_dash[m,i] = K_ads[i]/MM[i]*x[i]*gamma[i];
      C_forward[m,i] = max(1e-9,a_dash[m,i])^powerCoeffForward[m,i];
      C_reverse[m,i] = max(1e-9,a_dash[m,i])^powerCoeffReverse[m,i];
    end for;
    r_forward[m] = k_forward[m]*product(C_forward[m,:]);
    test_pCf[m] = product(C_forward[m,:]);
    r_reverse[m] = r_forward[m] - k_reverse[m]*product(C_reverse[m,:]);
    r_back[m] = k_reverse[m]*product(C_reverse[m,:]);
    test_pCr[m] = product(C_reverse[m,:]);
    sum_a_dash[m]=sum(a_dash[m,:]);
    r[m] =  m_cat/n* r_reverse[m]/(sum_a_dash[m]^2);
    K_bulk[m] = product(C_reverse[m,:])/product(C_forward[m,:]);
    PI[m].u = -K_bulk[m] + K_eq[m];
    for i in 1:nS loop
      if equilibrium[m] then
         NdotR[m,i] =0;
        else
 NdotR[m,i] =0*(1-omega)+omega*nu[m,i]/abs(nu[m,reacComp[m]])*r[m];
 end if;

  end for;
    end for;
  for i in 1:nS loop
    Ndot[i] = sum(NdotR[:,i]);
  end for;

/*** reaction enthalpy ***/

    for m in 1:nR loop
       h_R[m] =a_R[m]*T + b_R[m];
      //if reaction is exotherm, deltaH_R0 is negative, but we need an energy input in the energy equation
      // deltaH_singleReac[j,m] =- abs(NdotR[j,m,reacComp[m]])*deltaH_R0[m]/abs(nu[m,reacComp[m]]);
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
end ReactionBulkActivity;
