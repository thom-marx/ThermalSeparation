within ThermalSeparation.Media;
package Propane_Pentane_Liq "Propane, Pentane"

           constant Real phi[nSubstance]= {1,1}
  "association factor of each substance, if this substance is to be the solvent - used for claculation of diffusion coeffcients";


  extends ThermalSeparation.Media.BaseMediumLiquid(
                                             has_etaSubstance={true, true}, Tcrit= {
     369.8,
     469.6}, pcrit= {
     42.4e5,
     33.7e5},  Vcrit = {
    203e-6,
    304e-6}, omega = {
  0.153,
    0.251}, mu = {
  0,
0}, MMX= { 0.0441,  0.07215}, eq_Tsonopoulos = {1,  1},
     henry={false,false});

       constant Integer n = nSubstance "number of components";
     constant SI.MolarMass MM_substances[n]=MMX;


  redeclare replaceable model extends BaseProperties
      CalcSpecificEnthalpy calcSpecificEnthalpy(T0=T0, p=p, T=T, x=x);
  Real R;
protected
    SI.Temperature T_sat[nSubstance]={231.1,309.2};

    final parameter Integer counter[6]={0, 7, 14, 21, 28, 35};

    //dynamic viscosity for each component in the liquid phase
    //if a pure component does not exist as liquid at system pressure and temperature any value can be given
    //as long as the corresponding entry in the constant eta_pure is false

protected
    Real Theta = min(1,T/647.14) "dimensionless temperature";

  /* saturation pressure */

  Real A[nSubstance]={-6.71791,-7.27889};
  Real B[nSubstance]={1.33932,1.48754};
  Real C[nSubstance]={-2.23017,-2.87593};
  Real D[nSubstance]={-1.2499,-2.22146};

  equation
     for i in 1:nSubstance loop
     log(p_sat[i]/pcrit[i])=Tcrit[i]/T*(A[i]*(1-T/Tcrit[i])+B[i]*(1-T/Tcrit[i])^1.5+C[i]*(1-T/Tcrit[i])^3+D[i]*(1-T/Tcrit[i])^6);
     end for;

    //sigma =  235.8e-3*(1 - Theta)^1.256*(1 - 0.625*(1 - Theta));//surfaceTension(sat);
  sigma =  0.016;
    MM = max(1e-5, sum(x.*MM_substances));

          /*** density ***/
        d = 528*x[1]+646*x[2];

    /*** viscosity ***/

    eta = 0.00011*x[1]+0.00024*x[2];

    R = Modelica.Constants.R/MM;

     h = calcSpecificEnthalpy.h;
     u = h- p/d*MM;
       eta_comp= fill(eta,nSubstance);
    //d_water = 1000; //dummy
    v = MM/d;
      lambda = 0.108*x[1]+0.122*x[2];
      cp=x[1]*2456+x[2]*2206;

    annotation (structurallyIncomplete);
  end BaseProperties;


        redeclare model extends CalcSpecificEnthalpy
          /*** enthalpy ***/

        SI.Temperature T_sat[nSubstance]={231.1,309.2};
        equation
        /*** enthalpy ***/

          //  h = x[1]*2456*MMX[1]*(T_sat[1]-0)+x[2]*2206*MMX[2]*(T_sat[2]-0);
        // h = (x[1]*2456*MMX[1]+x[2]*2206*MMX[2])*(T-0);

        h = (x[1]*2456*MMX[1]+x[2]*2206*MMX[2])*(T-1392);

        //reference temperature T0=1392 has been modified so that h matches the value of FP

          annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
        end CalcSpecificEnthalpy;


redeclare model extends HenryCoefficient
     parameter SI.Pressure henry_H[nSubstance] = {1.63e8,0} "Henry coefficient";
     parameter Real henry_C[nSubstance]={2400,0}
    "constant to calculate temperature dependency";
     parameter SI.Temperature henry_T_norm=298 "norm temperature";
     parameter Boolean henry_temp = true;

      replaceable model HenryCoeff =
      ThermalSeparation.Media.Correlations.HenryCoefficient.Exponential                                   constrainedby
    ThermalSeparation.Media.Correlations.HenryCoefficient.BaseHenry;
    HenryCoeff henryCoeff(nS=nSubstance, henry_temp=henry_temp, henry_C=henry_C, henry_H=henry_H, x_l=x_l, T=T);
equation
  for i in 1:nSubstance loop
    if henry[i] then
    He[i] = henryCoeff.He[i];
    end if;
  end for;
end HenryCoefficient;


    redeclare model extends DiffusionCoefficient
    /*** model which shall be used to calculate the diffusion coefficients in a binary mixture ***/
      model D_Molecules =
      ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.Molecule.Wilke_Chang_aq
      (      phi=phi, MMX=MMX);

    /*** model where the diffusion coefficients in a binary mixture are used to calculate the diffusion coefficients in a multicomponent mixture ***/
    ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.DiffMolecules
      diffCoeff(
              redeclare model D_Molecules =
        D_Molecules,
    T=T, p=p, x=x, nS=nSubstance, has_etaSubstance=has_etaSubstance, eta=eta);
    equation
    D = diffCoeff.D;
    end DiffusionCoefficient;


  redeclare replaceable model extends ActivityCoefficient
  "calculates the activity coefficient for each component in the liquid phase"

    replaceable model ActivityCoeff =
        ThermalSeparation.Media.Correlations.ActivityCoefficient.Ideal                                      constrainedby
    ThermalSeparation.Media.Correlations.ActivityCoefficient.BaseActivityCoefficient
      annotation (choicesAllMatching=true, Dialog(tab="General", group=
            "Thermodynamic Equilibrium"));

    ActivityCoeff activityCoeff(nS=nSubstance, T=T, x_l = x_l);

  equation
     gamma = activityCoeff.gamma;

    annotation (structurallyIncomplete);
  end ActivityCoefficient;


  redeclare replaceable model extends FugacityCoefficient
  "calculates the fugacity coefficient at the saturation point for each component in the liquid phase"

      replaceable model FugacityCoeff =
      ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient.SaturationFugacitycoefficient
                                                                                                    constrainedby
    ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient.BaseFugacityCoefficient
      annotation (choicesAllMatching=true, Dialog(tab="General", group=
            "Thermodynamic Equilibrium"));

       FugacityCoeff fugacityCoeff(nS=nSubstance, T=T,p=p, p_sat=p_sat, Tcrit=Tcrit, pcrit=pcrit, omega=omega, NoOfEq=eq_Tsonopoulos, mu=mu);

  equation
    phi_sat = fugacityCoeff.phi_sat;
    annotation (structurallyIncomplete);
  end FugacityCoefficient;


  redeclare replaceable model extends ThermodynamicFactor
  equation
  Gamma=diagonal(ones(nSubstance));
    annotation (structurallyIncomplete);
  end ThermodynamicFactor;


  annotation (Icon(graphics={            Rectangle(
          extent={{-80,100},{100,-80}},
          lineColor={0,0,0},
          fillColor={215,230,240},
          fillPattern=FillPattern.Solid), Rectangle(
          extent={{-100,80},{80,-100}},
          lineColor={0,0,0},
          fillColor={240,240,240},
          fillPattern=FillPattern.Solid)}),
                          Documentation(info="<HTML>
<p>
This model calculates medium properties 
for water in the <b>liquid</b>, <b>gas</b> and <b>two phase</b> regions 
according to the IAPWS/IF97 standard, i.e., the accepted industrial standard
and best compromise between accuracy and computation time.
For more details see <a href=\"Modelica://Modelica.Media.Water.IF97_Utilities\">
Modelica.Media.Water.IF97_Utilities</a>. Three variable pairs can be the 
independent variables of the model:
<ol>
<li>Pressure <b>p</b> and specific enthalpy <b>h</b> are the most natural choice for general applications. This is the recommended choice for most general purpose applications, in particular for power plants.</li>
<li>Pressure <b>p</b> and temperature <b>T</b> are the most natural choice for applications where water is always in the same phase, both for liquid water and steam.</li>
<li>Density <b>d</b> and temperature <b>T</b> are explicit variables of the Helmholtz function in the near-critical region and can be the best choice for applications with super-critical or near-critial states.</li>
</ol>
The following quantities are always computed:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td valign=\"top\"><b>Variable</b></td>
      <td valign=\"top\"><b>Unit</b></td>
      <td valign=\"top\"><b>Description</b></td></tr>
  <tr><td valign=\"top\">T</td>
      <td valign=\"top\">K</td>
      <td valign=\"top\">temperature</td></tr>
  <tr><td valign=\"top\">u</td>
      <td valign=\"top\">J/kg</td>
      <td valign=\"top\">specific internal energy</b></td></tr>
  <tr><td valign=\"top\">d</td>
      <td valign=\"top\">kg/m^3</td>
      <td valign=\"top\">density</td></tr>
  <tr><td valign=\"top\">p</td>
      <td valign=\"top\">Pa</td>
      <td valign=\"top\">pressure</td></tr>
  <tr><td valign=\"top\">h</td>
      <td valign=\"top\">J/kg</td>
      <td valign=\"top\">specific enthalpy</b></td></tr>
</table>
<p>
In some cases additional medium properties are needed.
A component that needs these optional properties has to call
one of the functions listed in 
<a href=\"Modelica://Modelica.Media.UsersGuide.MediumUsage.OptionalProperties\">
Modelica.Media.UsersGuide.MediumUsage.OptionalProperties</a> and in
<a href=\"Modelica://Modelica.Media.UsersGuide.MediumUsage.TwoPhase\">
Modelica.Media.UsersGuide.MediumUsage.TwoPhase</a>.
</p>
</p>
<p>Many further properties can be computed. Using the well-known Bridgman's Tables, all first partial derivatives of the standard thermodynamic variables can be computed easily.</p>
</HTML>
"));
end Propane_Pentane_Liq;
