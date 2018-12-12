within ThermalSeparation.Media.WaterBasedLiquid;
package CO2_H2O "water: CO2, pure H2O for Siemens PilotPlant"
//adapted dynamic viscosity and heat capacity
constant Real f_psat = 1;
constant Boolean psat_Antoine = false;
       constant Real phi[nSubstance]= {1,2.26}
    "association factor of each substance, if this substance is to be the solvent - used for claculation of diffusion coeffcients";

  extends
    ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.PartialWaterBased(has_etaSubstance={ false, true}, Tcrit= {
     304.12,
      647.096}, pcrit= {
     7.374e6,
     2.2064e7}, Vcrit = {
  94.07e-6,
    55.95e-6}, omega = {
  0.225,
    0.344}, mu = { 0,
1.8},  MMX= { 0.0440095,  0.01801528}, eq_Tsonopoulos = { 1,  6},
     mediumName="WaterIF97",
     substanceNames={"water"},
     final reducedX=true,
     singleState=false,
     SpecificEnthalpy(start=1.0e5, nominal=5.0e5),
     Density(start=150, nominal=500),
     AbsolutePressure(start=50e5, nominal=10e5),
     Temperature(start=500, nominal=500),
     smoothModel=true,
     onePhase=true,
     fluidConstants=ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.waterConstants,
     nSubstance= 2,
     henry={true,false});

 // N2,
 // H2O,
 // O2,
 // CO2,
 // SO2,
 // HCl,
 // HF
 //
 //
 // Tcrit
 //    126.2,
 //    647.14,
 //    154.58,
 //    304.12,
 //    430.8,
 //    324.6,
 //    461.15
 //
 //  pcrit
 //    3.398e6,
 //    2.2064e7,
 //    5.043e6,
 //    7.374e6,
 //    7.884e6,
 //    8.31e6,
 //    6.48e6,
 //
 // Vcrit
 // 90.1e-6,
 // 55.95e-6,
 // 73.37e-6,
 // 94.07e-6,
 // 122e-6,
 // 222e-6,
 // 3448e-6,
 //
 // omega
 // 0.037,
 // 0.344,
 // 0.022,
 // 0.225,
 // 0.245,
 // 0.131,
 // 0.329
 //
 // mu
 // 0,
 // 1.8
 // 0,
 // 0,
 // 1.6
 // 0,
 // 0

       constant Integer n = nSubstance "number of components";
     constant SI.MolarMass MM_substances[n]=MMX;

  redeclare replaceable model extends BaseProperties(
    h(stateSelect=if ph_explicit and preferredMediumStates then StateSelect.prefer else StateSelect.default),
    d(stateSelect=if dT_explicit and preferredMediumStates then StateSelect.prefer else StateSelect.default),
    T(stateSelect=if (pT_explicit or dT_explicit) and preferredMediumStates then StateSelect.prefer else StateSelect.default),
    p(stateSelect=if (pT_explicit or ph_explicit) and preferredMediumStates then StateSelect.prefer else StateSelect.default))
    "Base properties of water"

    SaturationProperties sat(Tsat(start=300.0), psat(start=1.0e5))
      "saturation temperature and pressure";

          CalcSpecificEnthalpy calcSpecificEnthalpy(T0=T0, p=p, T=T, x=x);

  /***density and enthalpy ***/
  protected
    Real pi1 "dimensionless pressure";
    Real tau1 "dimensionless temperature";
    Real[45] o "vector of auxiliary variables";
    Real pi;
    Real tau;
    Real gtau;

  /*** viscosity ***/
    constant Real n0=1.0 "viscosity coefficient";
    constant Real n1=0.978197 "viscosity coefficient";
    constant Real n2=0.579829 "viscosity coefficient";
    constant Real n3=-0.202354 "viscosity coefficient";
    constant Real[42] nn=array(0.5132047,0.3205656, 0.0, 0.0, -0.7782567,
        0.1885447, 0.2151778, 0.7317883, 1.241044, 1.476783, 0.0, 0.0, -0.2818107,
         -1.070786, -1.263184, 0.0, 0.0, 0.0, 0.1778064, 0.460504,
        0.2340379, -0.4924179, 0.0, 0.0, -0.0417661, 0.0, 0.0, 0.1600435,
        0.0, 0.0, 0.0, -0.01578386, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.003629481,
         0.0, 0.0) "viscosity coefficients";
    constant SI.Density rhostar=317.763 "scaling density";
    constant SI.DynamicViscosity etastar=55.071e-6 "scaling viscosity";
    constant SI.Temperature tstar=647.226 "scaling temperature";
    Real delta "dimensionless density";
    Real deltam1 "dimensionless density";
    Real tau_eta "dimensionless temperature";
    Real taum1 "dimensionless temperature";
    Real Psi0 "auxiliary variable";
    Real Psi1 "auxiliary variable";
    Real tfun[6] "auxiliary variable";
    Real rhofun[7] "auxiliary variable";
    Real test[42];
    final parameter Integer counter[6]={0, 7, 14, 21, 28, 35};

    //dynamic viscosity for each component in the liquid phase
    //if a pure component does not exist as liquid at system pressure and temperature any value can be given
    //as long as the corresponding entry in the constant eta_pure is false

  public
    SI.SpecificHeatCapacity cp_solution= 2.8571*T + 2505.7;

  /*** surface tension ***/
  protected
    Real Theta = min(1,T/647.14) "dimensionless temperature";
  /*** saturation pressure p_sat ***/
    Real[8] o_s "vector of auxiliary variables";
    Real Tlim=min(T, 647.096);
    Real test2=max(0,if time< 1 then 1 else 1-(time-1)/100);

  equation
    sigma =  235.8e-3*(1 - Theta)^1.256*(1 - 0.625*(1 - Theta));//surfaceTension(sat);
  MM=max(1e-5, sum(x.*MM_substances));

  /*** enthalpy ***/
        h=calcSpecificEnthalpy.h;//x[2]*4180*MM_substances[2]*(T-0) + x[1]*21e3;

          /*** density ***/
        d = d_water;
      pi = p/data.PSTAR1;
      pi1 = 7.1000000000000 - pi;
      tau = data.TSTAR1 /T;
      tau1 = -1.22200000000000 + tau;
        o[1] =  tau1*tau1;
    o[2] = o[1]*o[1];
    o[3] = o[2]*o[2];
    o[4] = o[3]*tau1;
    o[5] = 1/o[4];
    o[6] = o[1]*o[2];
    o[7] = o[1]*tau1;
    o[8] = 1/o[7];
    o[9] = o[1]*o[2]*o[3];
    o[10] = 1/o[2];
    o[11] = o[2]*tau1;
    o[12] = 1/o[11];
    o[13] = o[2]*o[3];
    o[14] = 1/o[3];
    o[15] = pi1*pi1;
    o[16] = o[15]*pi1;
    o[17] = o[15]*o[15];
    o[18] = o[17]*o[17];
    o[19] = o[17]*o[18]*pi1;
    o[20] = o[15]*o[17];
    o[21] = o[3]*o[3];
    o[22] = o[21]*o[21];
    o[23] = o[22]*o[3]*tau1;
    o[24] = 1/o[23];
    o[25] = o[22]*o[3];
    o[26] = 1/o[25];
    o[27] = o[1]*o[2]*o[22]*tau1;
    o[28] = 1/o[27];
    o[29] = o[1]*o[2]*o[22];
    o[30] = 1/o[29];
    o[31] = o[1]*o[2]*o[21]*o[3]*tau1;
    o[32] = 1/o[31];
    o[33] = o[2]*o[21]*o[3]*tau1;
    o[34] = 1/o[33];
    o[35] = o[1]*o[3]*tau1;
    o[36] = 1/o[35];
    o[37] = o[1]*o[3];
    o[38] = 1/o[37];
    o[39] = 1/o[6];
    o[40] = o[1]*o[22]*o[3];
    o[41] = 1/o[40];
    o[42] = 1/o[22];
    o[43] = o[1]*o[2]*o[21]*o[3];
    o[44] = 1/o[43];
    o[45] = 1/o[13];

      gtau =  pi1*(o[38]*(-0.00254871721114236 + o[1]*(0.0042494411096112
       + (0.0189900682184190 + (-0.0218417171754140 - 0.000158515073909790*
      o[1])*o[1])*o[6])) + pi1*(o[10]*(0.00141552963219801 + o[2]*(
      0.000047661393906987 + o[1]*(-0.0000132425535992538 -
      1.23581493705910e-14*o[9]))) + pi1*(o[12]*(0.000126718579380216 -
      5.1123076872062e-9*o[37]) + pi1*(o[39]*(0.0000112126409540000 + (
      1.30342445791202e-6 - 1.43417299379240e-12*o[13])*o[7]) + pi1*(
      3.2413597488094e-6*o[5] + o[16]*((1.40077319158051e-8 +
      1.04549227383804e-9*o[11])*o[45] + o[19]*(1.99410180757040e-17*o[44]
       + o[15]*(-4.4882754268415e-19*o[42] + o[20]*(-1.00075970318621e-21*o[
      28] + pi1*(4.6595728296277e-22*o[26] + pi1*(-7.2912378325616e-23*o[24]
       + 3.8350205789908e-24*o[41]*pi1))))))))))) + o[8]*(-0.292659424263340
       + tau1*(0.84548187169114 + o[1]*(3.3855169168385 + tau1*(-1.91583926775744
       + tau1*(0.47316115539684 + (-0.066465668798004 + 0.0040607314991784*
      tau1)*tau1)))));

      sat.psat = p;
      sat.Tsat = saturationTemperature(p);

     /*** viscosity ***/
     delta = d/rhostar;
    deltam1 = delta - 1.0;
    tau_eta = tstar/T;
    taum1 =  tau_eta - 1.0;
    Psi0 =  1/(n0 + (n1 + (n2 + n3*tau_eta)*tau_eta)*tau_eta)/(tau_eta^0.5);
    tfun[1] = 1;
    for k in 2:6 loop
      tfun[k] =  taum1^(k-1);
    end for;
    rhofun[1]= 1;
    for k in 2:7 loop
      rhofun[k] = deltam1^(k - 1);
    end for;
    for m in 1:6 loop
      for n in 0:6 loop
        test[counter[m]+n+1] = nn[m + n*6]*tfun[m]*rhofun[n + 1];
      end for;
    end for;
    Psi1 = sum(test);
    eta = etastar*Psi0*Modelica.Math.exp(delta*Psi1);

    /*** saturation pressure p_sat of water ***/
    o_s[1] = -650.17534844798 + Tlim;
    o_s[2] = 1/o_s[1];
    o_s[3] = -0.238555575678490*o_s[2];
    o_s[4] = o_s[3] + Tlim;
    o_s[5] = -4823.2657361591*o_s[4];
    o_s[6] = o_s[4]*o_s[4];
    o_s[7] = 14.9151086135300*o_s[6];
    o_s[8] = 405113.40542057 + o_s[5] + o_s[7];

       p_sat = fill(1,n).*fill(16.0e6*o_s[8]*o_s[8]*o_s[8]*o_s[8]*1/(3.2325550322333e6 -
         12020.8247024700*o_s[4] + 17.0738469400920*o_s[6] + (-4.0*(-724213.16703206
          + 1167.05214527670*o_s[4] + o_s[6])*o_s[8] + (-3.2325550322333e6 +
         12020.8247024700*o_s[4] - 17.0738469400920*o_s[6])^2.0)^0.5)^4.0,n);

    u = h - p/d*MM;
    R = Modelica.Constants.R/fluidConstants[1].molarMass;

    annotation (structurallyIncomplete);
  end BaseProperties;

  redeclare model extends CalcSpecificEnthalpy
    /*** enthalpy ***/
    parameter SI.SpecificHeatCapacity cp_water=4200;

  equation
  /*** enthalpy ***/
    h = cp_water*sum(x.*MM_substances)*(T-T0);
    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;

  model EvaporationEnthalpy
  /*** enthalpy ***/
  input SI.Pressure p;
  input SI.Temperature T;
  //input SI.MoleFraction x[n];
  ThermalSeparation.Units.MolarEnthalpy h[n];

  protected
  SI.SpecificEnthalpy r_water = -2462.5 * T +3177.8e3;
  equation
        /*** enthalpy ***/

  h[1] = 21e3;
  h[2] = MM_substances[2]*r_water;
    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end EvaporationEnthalpy;

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
        (    phi=phi, MMX=MMX);

    /*** model where the diffusion coefficients in a binary mixture are used to calculate the diffusion coefficients in a multicomponent mixture ***/
    ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.DiffMolecules
      diffCoeff(
              redeclare model D_Molecules =
        D_Molecules,
    T=T, p=p, x=x, nS=nSubstance, has_etaSubstance=has_etaSubstance, eta=eta);
    equation
    D = diffCoeff.D;
    end DiffusionCoefficient;

redeclare model extends SaturationPressure

equation
  /*** saturation pressure p_sat ***/
  o_s[1] = -650.17534844798 + Tlim;
  o_s[2] = 1/o_s[1];
  o_s[3] = -0.238555575678490*o_s[2];
  o_s[4] = o_s[3] + Tlim;
  o_s[5] = -4823.2657361591*o_s[4];
  o_s[6] = o_s[4]*o_s[4];
  o_s[7] = 14.9151086135300*o_s[6];
  o_s[8] = 405113.40542057 + o_s[5] + o_s[7];

// if psat_Antoine then
//    /*** Antoine Siemens ***/
//  p_sat = fill(exp(24.9022-5342.9/(T+21.2826)),n);
//  else
      p_sat = fill(f_psat,n).*fill(16.0e6*o_s[8]*o_s[8]*o_s[8]*o_s[8]*1/(3.2325550322333e6 -
        12020.8247024700*o_s[4] + 17.0738469400920*o_s[6] + (-4.0*(-724213.16703206
         + 1167.05214527670*o_s[4] + o_s[6])*o_s[8] + (-3.2325550322333e6 +
        12020.8247024700*o_s[4] - 17.0738469400920*o_s[6])^2.0)^0.5)^4.0,n);

  //   p_sat = fill(f_psat,n);
//end if;
end SaturationPressure;
  annotation (Icon(graphics={Text(
            extent={{-94,84},{94,40}},
            lineColor={127,191,255},
            textString=
               "IF97"), Text(
            extent={{-94,20},{94,-24}},
            lineColor={127,191,255},
            textString=
               "water")}),Documentation(info="<HTML>
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
end CO2_H2O;
