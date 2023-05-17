within ThermalSeparation.FilmModel.BaseClasses;
package EnhancementFactor
  "This package contains correlations to calculate the enhancement factor"
model IrrFirstOrder "enhancement factor: irreversible first order reaction"
  /*** A -> C ***/

  extends ThermalSeparation.FilmModel.BaseClasses.EnhancementFactor.BaseEnhancement;

  parameter Real a_R = 38.814 "h_R [J/mol] = a_R * T [K] + b_R";
  parameter Real b_R = -89155 "h_R [J/mol] = a_R * T [K] + b_R";

//reacComp corresponds to component A in van Swaaij and Versteeg
 parameter Integer reacComp = 1
      "index of component for which enhancement factor is calculated";
    parameter Boolean fast=false
      "reaction is known to be fast compared to mass transport";

  MediumLiquid.ReactionRates reactionRates(c=c, gamma=gamma, propsLiq=propsLiq);

  Real Ha "Hatta number";

MediumLiquid.DiffusionCoefficient diffCoeff(T=T, p=p, x=x, eta=eta_comp);
  protected
  SI.DiffusionCoefficient D[nS] = diffCoeff.D_diluted;

equation
  for i in 1:nS loop
    assert(D[i] > 0, "The medium model must supply values for the diluted diffusion coefficient, D_diluted!");
  end for;

    //equation 17 in van Swaaij and Versteeg: Mass Transfer Accompanied Wich Complex Reversible Chemical Reactions
    Ha =sqrt(reactionRates.k_forward[1]*D[reacComp])/k_l[reacComp];
    for i in 1:nS loop
      if i==reacComp then
        E[i]=if fast == true then Ha else Ha/Modelica.Math.tanh(Ha);
      else
        E[i]=1;
      end if;
      end for;

    annotation (Documentation(info="<html>
<p>reaction type: <code><font style=\"color: #006400; \">A&nbsp;-&GT;&nbsp;C</font></code></p>
<pre><font style=\"color: #006400; \">Enhancement factor is calculated for component A.</font></pre>
</html>"));
end IrrFirstOrder;

model IrrPseudoFirstOrder "enhancement factor: irreversible pseudo first order"
  /*** A + B -> C + D mit c_A << c_B ***/
  extends ThermalSeparation.FilmModel.BaseClasses.EnhancementFactor.BaseEnhancement;

  parameter Real a_R = 38.814 "h_R [J/mol] = a_R * T [K] + b_R";
  parameter Real b_R = -89155 "h_R [J/mol] = a_R * T [K] + b_R";

//reacComp corresponds to component A in van Swaaij and Versteeg
 parameter Integer reacComp = 1
      "index of component for which enhancement factor is calculated";
    parameter Integer reacCompB = 2 "index of component reacting with A";
 constant Integer nR=MediumLiquid.nR "number of reactions";
   constant Real nu[nR,nS]=MediumLiquid.nu
      "stocheometric coefficients for the reaction";

   MediumLiquid.ReactionRates reactionRates(c=c, gamma=gamma, propsLiq=propsLiq);

  Real Ha "Hatta number";
  Real E_max;

MediumLiquid.DiffusionCoefficient diffCoeff(T=T, p=p, x=x, eta=eta_comp);
  protected
  SI.DiffusionCoefficient D[nS] = diffCoeff.D_diluted;

equation
  for i in 1:nS loop
    assert(D[i] > 0, "The medium model must supply values for the diluted diffusion coefficient, D_diluted!");
  end for;

    //equation 31 in van Swaaij and Versteeg: Mass Transfer Accompanied Wich Complex Reversible Chemical Reactions; with m=1 and n=1
    Ha =sqrt(reactionRates.k_forward[1]*c[reacCompB]*D[reacComp])/k_l[reacComp];
    E_max = 1+D[reacCompB]*c[reacCompB]/(abs(nu[1,reacCompB])*D[reacComp]*c_l_star[reacComp]);
    for i in 1:nS loop
      if i==reacComp then
        E[i]= Ha*sqrt((E_max-E[reacComp])/(E_max-1))/(tanh(Ha*sqrt((E_max-E[reacComp])/(E_max-1))));
      else
        E[i]=1;
      end if;
      end for;

    annotation (Documentation(info="<html>
<p><code>reaction type: </code><code><font style=\"color: #006400; \">A&nbsp;+&nbsp;B&nbsp;-&GT;&nbsp;C&nbsp;+&nbsp;D&nbsp;with&nbsp;c_A&nbsp;&LT;&LT;&nbsp;c_B</font></code></p>
<font style=\"color: #006400; \">Enhancement factor is calculated for component A.</font></pre>
</html>"));
end IrrPseudoFirstOrder;

model RevAtoC "enhancement factor: reversible reaction,  A <-> C"
  /*** A <-> C ***/
  //film model, instanteneous reaction, only diluted systems (Olander 1960)
  extends ThermalSeparation.FilmModel.BaseClasses.EnhancementFactor.BaseEnhancement;

  parameter Real a_R = 38.814 "h_R [J/mol] = a_R * T [K] + b_R";
  parameter Real b_R = -89155 "h_R [J/mol] = a_R * T [K] + b_R";

//reacComp corresponds to component A in van Swaaij and Versteeg
 parameter Integer reacComp = 1
      "index of component for which enhancement factor is calculated";
    parameter Integer reacCompC = 2
      "index of component which is the reaction product";
  constant Integer nR=MediumLiquid.nR "number of reactions";
  constant Real nu[nR,nS]=MediumLiquid.nu
      "stocheometric coefficients for the reaction";
   MediumLiquid.EquilibriumConstant eqConst(c=c, gamma=gamma, propsLiq=propsLiq);

MediumLiquid.DiffusionCoefficient diffCoeff(T=T, p=p, x=x, eta=eta_comp);
  protected
  SI.DiffusionCoefficient D[nS] = diffCoeff.D_diluted;

equation
  for i in 1:nS loop
    assert(D[i] > 0, "The medium model must supply values for the diluted diffusion coefficient, D_diluted!");
      assert(nR==1, "this enhancement model is only valid if there is only one reaction in the liquid phase");
  end for;

 //equation can be found in van Swaaij and Versteeg (eq. 38) or Chang and Rochelle (table 1)
    for i in 1:nS loop
      if i==reacComp then
        E[i]= 1 + D[reacCompC]/D[reacComp]*eqConst.K_eq[1];
      else
        E[i]=1;
      end if;
      end for;

    annotation (Documentation(info="<html>
<p><code>reaction type: </code><code><font style=\"color: #006400; \">A&nbsp;&LT;-&GT;&nbsp;C</font></code></p>
<font style=\"color: #006400; \">Enhancement factor is calculated for component A.</font></pre>
</html>"));
end RevAtoC;

model RevAto2C "enhancement factor: reversible reaction,  A <-> 2 C"
  /*** A <-> 2C ***/
  //film model, instanteneous reaction, only diluted systems (Olander 1960)
  extends ThermalSeparation.FilmModel.BaseClasses.EnhancementFactor.BaseEnhancement;

    parameter Real a_R = 38.814 "h_R [J/mol] = a_R * T [K] + b_R";
  parameter Real b_R = -89155 "h_R [J/mol] = a_R * T [K] + b_R";

//reacComp corresponds to component A in van Swaaij and Versteeg
 parameter Integer reacComp = 1
      "index of component for which enhancement factor is calculated";
    parameter Integer reacCompC = 2
      "index of component which is the reaction product";

   constant Integer nR=MediumLiquid.nR "number of reactions";
  constant Real nu[nR,nS]=MediumLiquid.nu
      "stocheometric coefficients for the reaction";
   MediumLiquid.EquilibriumConstant eqConst(c=c, gamma=gamma, propsLiq=propsLiq);

MediumLiquid.DiffusionCoefficient diffCoeff(T=T, p=p, x=x, eta=eta_comp);
  protected
  SI.DiffusionCoefficient D[nS] = diffCoeff.D_diluted;

equation
  for i in 1:nS loop
    assert(D[i] > 0, "The medium model must supply values for the diluted diffusion coefficient, D_diluted!");
    assert(nR==1, "this enhancement model is only valid if there is only one reaction in the liquid phase");
  end for;

 //equation can be found in van Swaaij and Versteeg (eq. 39) or Chang and Rochelle (table 1)
    for i in 1:nS loop
      if i==reacComp then
        E[i]= 1+D[reacCompC]/D[reacComp]*sqrt(0.25*eqConst.K_eq[1])/(sqrt(c_l_star[reacComp])+sqrt(c[reacComp]));
      else
        E[i]=1;
      end if;
      end for;

    annotation (Documentation(info="<html>
<p><code>reaction type: </code><code><font style=\"color: #006400; \">&nbsp;A&nbsp;&LT;-&GT;&nbsp;2&nbsp;C</font></code></p>
<font style=\"color: #006400; \">Enhancement factor is calculated for component A.</font></pre>
</html>"));
end RevAto2C;

model RevAplusBtoC "enhancement factor: reversible reaction,  A + B <-> C"
/*** A + B <-> C ***/
  //film model, instanteneous reaction, only diluted systems (Olander 1960)
  extends ThermalSeparation.FilmModel.BaseClasses.EnhancementFactor.BaseEnhancement;
    parameter Real a_R = 38.814 "h_R [J/mol] = a_R * T [K] + b_R";
  parameter Real b_R = -89155 "h_R [J/mol] = a_R * T [K] + b_R";

    //reacComp corresponds to component A in van Swaaij and Versteeg
 parameter Integer reacComp = 1
      "index of component for which enhancement factor is calculated";
        parameter Integer reacCompB = 2 "index of component reacting with A";
    parameter Integer reacCompC = 2
      "index of component which is the reaction product";
   constant Integer nR=MediumLiquid.nR "number of reactions";
  constant Real nu[nR,nS]=MediumLiquid.nu
      "stocheometric coefficients for the reaction";
   MediumLiquid.EquilibriumConstant eqConst(c=c, gamma=gamma, propsLiq=propsLiq);

MediumLiquid.DiffusionCoefficient diffCoeff(T=T, p=p, x=x, eta=eta_comp);
  protected
  SI.DiffusionCoefficient D[nS] = diffCoeff.D_diluted;

equation
  for i in 1:nS loop
    assert(D[i] > 0, "The medium model must supply values for the diluted diffusion coefficient, D_diluted!");
      assert(nR==1, "this enhancement model is only valid if there is only one reaction in the liquid phase");
  end for;

 //equation can be found in van Swaaij and Versteeg (eq. 40) or Chang and Rochelle (table 1)
    for i in 1:nS loop
      if i==reacComp then
        E[i]= 1+D[reacCompB]*c[reacCompB]*eqConst.K_eq[1]/(D[reacComp]*(c_l_star[reacComp]*eqConst.K_eq[1] + D[reacCompB]/(D[reacCompC])));
      else
        E[i]=1;
      end if;
      end for;

    annotation (Documentation(info="<html>
<pre><font style=\"color: #006400; \">Reversible&nbsp;reaction,&nbsp;&nbsp;A&nbsp;+&nbsp;B&nbsp;&LT;-&GT;&nbsp;C</font>
<font style=\"color: #006400; \">Enhancement factor is calculated for component A.</font></pre>
</html>"));
end RevAplusBtoC;

  partial model BaseEnhancement "partial model for enhancement factor"

    parameter Integer nS=2
      "number of components, educts=negative, products=positive" annotation(Dialog(enable=false));
    input SI.Temperature T;
   input SI.Concentration c[nS];
    input SI.Concentration c_l_star[nS];
    input SI.MoleFraction x[nS];
    input SI.Pressure p;
    input SI.DynamicViscosity eta_comp[nS];
    input ThermalSeparation.Units.CoefficentOfMassTransfer k_l[nS]
      "liquid side mass transfer coefficient";
    input Real gamma[nS];
    input MediumLiquid.ThermodynamicProperties propsLiq;

    replaceable package MediumLiquid = Media.BaseMediumLiquidReaction annotation(Dialog(enable=false));

    output Real E[ nS] "enhancement factor for each substance";

    annotation (Icon(graphics));
  end BaseEnhancement;

  model Constant "constant enhancement factor"
    extends BaseEnhancement;
     parameter Integer reacComp
      "index of component for which enhancement factor is calculated";
      parameter Real E_const "value for enhancement factor";
  equation
  for i in 1:nS loop
    E[i]=if i==reacComp then E_const else 1;
  end for;
  end Constant;
  annotation (Documentation(info="<html>
<p>This package contains some models to calculate the enhancement factor. Enhancement factors for parallel reactions are not implemented since it is not clear how to calculate the reaction enthalpy in this case.</p>
</html>"));
end EnhancementFactor;
