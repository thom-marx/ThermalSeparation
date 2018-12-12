within ThermalSeparation.Media;
package Methylacetatsynthese_Liq "Methylacetat, Essigsäure, Methanol, Wasser"
  extends BaseMediumLiquidReaction(
  has_etaSubstance=fill(true,nSubstance),
  nSubstance=4,
  Tcrit={506.1, 592.7, 513.1, 647.3},
  pcrit={4.6e6, 57.9e5, 8.084e6, 2.2064e7},
  Vcrit={228e-6, 171e-6, 118e-6, 55.95e-6},
  omega={ 0.326,0.447,0.565,0.022},
  mu={1.7,1.3,1.7, 1.8},
  MMX= {0.0741,0.0601,0.03204,0.018},
  eq_Tsonopoulos = {4,4,4, 6},
  henry=fill(false,nSubstance),
 nR=1, nu={{1,-1,-1,1}},  reacComp={1}, nS_reac=4, redeclare model
      MolarReactionEnthalpy =
        ThermalSeparation.Media.Correlations.Reaction.MolarReactionEnthalpy.EnthalpyOfFormation
        ( nR=nR));

  constant Integer n = nSubstance "number of components";

    constant SI.MolarMass MM_substances[n] = MMX;

  redeclare replaceable model extends BaseProperties(T0 = 298.15)

   SI.MassFraction w[n](start=fill(0.03,n));

  constant SI.DynamicViscosity sigma_EtOH=0.02255;// bei 20 C N/m
  constant SI.DynamicViscosity sigma_H2O=0.0728;// bei 20 C N/m

  ThermodynamicState state( T=T, h=h/MM, d=d, p=p)
      "thermodynamic state variables for optional functions";
  Real MM_min;
  Real MM_var;

     CalcSpecificEnthalpy calcSpecificEnthalpy(nS=n,T0=T0, p=p, T=T, x=x);

  /***density and enthalpy ***/
  DensityWater densityWater(p=1e5,T=T,scale=1);
  DensityWater densityMethanol(p=p,T=T,scale=0.8);

    Real d_Methanol=d3.rho;//densityMethanol.d;
      Real d_Essig=d2.rho;//1050;
      SI.Density d_Water=densityWater.d;
      SI.Density d_Methylacetate = d1.rho;//930;

      SaturationDensityDIPPR d1(T=T, C={1130,0.259,506.5,0.2764}, MM=MM_substances[1]);
          SaturationDensityDIPPR d2(T=T, C={1448.6, 0.2589, 591.95, 0.2529}, MM=MM_substances[2]);
          SaturationDensityInfochem d3(T=T, C={8474.58, 22277.5, 0.375}, MM=MM_substances[3], Tc=Tcrit[3]);

      Real d_molar_Methanol=d3_molar.rho;//Einheit: mol/m3;
      Real d_molar_Essig=d2_molar.rho;//1050;
      SI.Density d_molar_Water=densityWater.d/MM_substances[4];
      SI.Density d_molar_Methylacetate = d1_molar.rho;//930;
          MolarSaturationDensityDIPPR d1_molar(T=T, C={1130,0.259,506.5,0.2764}, MM=MM_substances[1]);
          MolarSaturationDensityDIPPR d2_molar(T=T, C={1448.6, 0.2589, 591.95, 0.2529}, MM=MM_substances[2]);
          MolarSaturationDensityInfochem d3_molar(T=T, C={8474.58, 22277.5, 0.375}, MM=MM_substances[3], Tc=Tcrit[3]);

  /*** saturation pressure ***/
    Real A_Ant[:] = {16.2685, 17.4068, 18.6073, 18.585};
    Real B_Ant[:] = {2665.5416, 3785.565, 3643.3113, 3984.9228};
    Real C_Ant[:] = {-53.424, -39.626, -33.4240, -39.724};

   Real shift = max(0,min((time-100)/100,0.9));

  equation
      eta_comp= fill(eta,n);
          for j in 1:n loop
    w[j] =MMX[j]/MM*x[j];
    end for;

    eta=0.001002;
    sigma=sigma_EtOH*x[1]+sigma_H2O*x[2];

  /*** density - ok ***/
  v=x[1]/d_molar_Methylacetate + x[2]/d_molar_Essig + x[3]/d_molar_Methanol + x[4]/d_molar_Water;
   d=1/((1/max(1e-5,d_Methanol))*max(1e-5,w[3]) + (1/max(1e-5,d_Essig))*max(1e-5,w[2]) + (1/max(1e-5,d_Methylacetate))*max(1e-5,w[1]) + (1/max(1e-5,d_Water))*max(1e-5,w[4]));
    MM_min=min(MMX);
    MM= MM_var;
    MM_var = max(MM_min,sum(x[j]*MMX[j] for j in 1:n));

  /*** specific enthalpy and inner energy -ok ***/
    h = calcSpecificEnthalpy.h;
   // u = h - p/d*MM;
   u = h - p*v;

  /*** saturation pressure -ok ***/
  for i in 1:n loop
    p_sat[i] = 133.3*exp(A_Ant[i] - B_Ant[i]/(T+C_Ant[i]));
  end for;

      lambda =0.55;
      cp=4000;

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end BaseProperties;

redeclare replaceable model extends ActivityCoefficient
  /*** activity coefficient -ok ***/
   parameter Real r_UNI[:] = {2.8042, 2.2024, 1.4311, 0.92};
   parameter Real q_UNI[:] = {2.576, 2.072, 1.432, 1.4};
//   parameter Real r_UNI[:] = {2.8042,1.9, 1.4311, 0.92};
//   parameter Real q_UNI[:] = {2.576, 1.8, 1.432, 1.4};
  Real a_UNI[:,:] = {{1, 195.08,  620.11, 593.7},    {-76.614, 1,390.26, 422.38},        {  -86.02, 65.245, 1, -235.52},  {    -265.83, -98.12, -761.48, 1}};
  Real b_UNI[:,:] = {{1, 0.40148,-0.96715, 0.010143},{ -0.33388, 1, 0.97039, -0.051007}, { 0.13943, -2.0346, 1, 0.41274}, { 0.96295, -0.29355, 5.2418, 1}};
  Real c_UNI[:,:] = {{1, 0, 0, -2.1609e-3},          {  0, 1, -3.0613e-3, -2.4019e-4},   { 0, 3.157e-3, 1, -1.5415e-3},   {   2.0113e-4, -7.6741e-5, -4.3663e-3, 1}};
  Real u_UNI[Methylacetatsynthese_Liq.n,Methylacetatsynthese_Liq.n];

 ThermalSeparation.Media.Correlations.ActivityCoefficient.UNIQUAC
      activityCoefficient(
                        nS=n,T=T,x_l=x_l,q=q_UNI, r=r_UNI, delta_u=u_UNI);
// Real omega = 0.5*tanh(0.01*(time-200))+0.5;
 Real omega = 0.5*tanh(0.08*(time-30))+0.5;
equation
  /*** activity coefficient -ok ***/
for i in 1:n loop
  for j in 1:n loop
    if i == j then
      u_UNI[i,j]=1;
    else
   //   if i==1 or j==1 then
        u_UNI[i,j] = (a_UNI[i,j] + b_UNI[i,j]*T + c_UNI[i,j]*T^2)*Modelica.Constants.R;
   //     end if;
    end if;
  end for;
      end for;
//        u_UNI[2,3]=3135;
//        u_UNI[3,2]=-2181;
//      u_UNI[2,4]=41.7*Modelica.Constants.R;
//      u_UNI[4,2]=-0.5244*Modelica.Constants.R;
//      u_UNI[3,4]=-165.3*Modelica.Constants.R;
//      u_UNI[4,3]=254.7*Modelica.Constants.R;

  for j in 1:n loop
    gamma[j] = 1*(1-omega) + omega*activityCoefficient.gamma[j];
  end for;

end ActivityCoefficient;

redeclare model extends FugacityCoefficient
    "calculates the fugacity coefficient at the saturation point for each component in the liquid phase"

  parameter Integer nS=4;
  parameter Real A = -17.374;
  parameter Real B = -7290;
  Real f0 = 1e5 "standard fugacity";
  Real K_D_star[nS] = exp(A-B/T)*p_sat/f0;
  Real sq[nS];
equation
  for i in 1:nS loop
    sq[i]= sqrt(1+4*K_D_star[i]);
  end for;

  phi_sat[2] = (sq[2]-1)/(2*K_D_star[2]);
     for i in 1:1 loop
    phi_sat[i] =  1;
     end for;
   for i in 3:4 loop
    phi_sat[i] =  1;
     end for;

//     replaceable model FugacityCoeff =
//       ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient.Ideal                    constrainedby
//     ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient.BaseFugacityCoefficient
//     annotation (choicesAllMatching=true, Dialog(tab="General", group=
//           "Thermodynamic Equilibrium"));
//
//      FugacityCoeff fugacityCoeff(nS=nSubstance, T=T,p=p, p_sat=p_sat, Tcrit=Tcrit, pcrit=pcrit, omega=omega, NoOfEq=eq_Tsonopoulos, mu=mu);
//
// equation
//   phi_sat = fugacityCoeff.phi_sat;

end FugacityCoefficient;

  model Testmodell
    parameter Integer n=1;
    parameter Integer nS=4;
  parameter SI.Temperature T = 340;
  parameter SI.Pressure p = 100000;
  parameter SI.MoleFraction x[nS] = {  0.469120,  0.345584,  7.297487e-2,  0.112321};//,0.05}};

  BaseProperties BP(
      T=T,
      p=p,
      x=x);
  end Testmodell;

  redeclare model extends CalcSpecificEnthalpy
  /*** enthalpy ***/
  parameter Integer nS=4;

  /*** Wärmekapazitäten durch Polynome 4. Ordnung austauschen ***/
  SI.MolarHeatCapacity cp_MeOAc = 150;
  SI.MolarHeatCapacity cp_HOAc = 650;
  SI.MolarHeatCapacity cp_MeOH = 120;
  SI.MolarHeatCapacity cp_H2O = 4187*MMX[4];
  SI.MolarHeatCapacity cp[nS] = {cp_MeOAc, cp_HOAc, cp_MeOH, cp_H2O};

  /*** Standardbildungsenthalpien ***/
  ThermalSeparation.Units.MolarEnthalpy h_f[nS] = {-442790, -484089, -238572, -285830};
  ThermalSeparation.Units.MolarEnthalpy h_substance[nS];

  equation
        /*** enthalpy ***/
  h_substance =h_f+  cp*(T-T0);
  h = sum(x.*h_substance);
    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;

  model DensityWater
    input SI.Pressure p;
    parameter Real scale=1;
    input SI.Temperature T;
    output SI.Density d;
  protected
    Real pi1 "dimensionless pressure";
    Real tau1 "dimensionless temperature";
    Real[45] o "vector of auxiliary variables";
    Real pi;
    Real gpi;
    Real tau;

      constant SI.SpecificHeatCapacity RH2O=461.526
      "specific gas constant of water vapour";
    constant SI.Temperature TSTAR1=1386.0
      "normalization temperature for region 1 IF97";
    constant SI.Pressure PSTAR1=16.53e6
      "normalization pressure for region 1 IF97";

  equation
            /*** density ***/

      d =scale* p/( RH2O * T * pi * gpi);
      pi = p/ PSTAR1;
      gpi = pi1*(pi1*(o[10]*(0.000095038934535162 + o[2]*(
        8.4812393955936e-6 + 2.55615384360309e-9*o[6])) + pi1*(o[12]*(
        8.9701127632000e-6 + (2.60684891582404e-6 + 5.7366919751696e-13*o[13])
        *o[7]) + pi1*(2.02584984300585e-6*o[14] + o[16]*((1.01874413933128e-8
         + 1.39398969845072e-9*o[11])*o[36] + o[19]*(1.44400475720615e-17*o[
        34] + o[15]*(-3.3300108005598e-19*o[32] + o[20]*(-7.6373766822106e-22
        *o[30] + pi1*(3.5842867920213e-22*o[28] + pi1*(-5.6507093202352e-23*o[
        26] + 2.99318679335866e-24*o[24]*pi1))))))))) + o[8]*(
        0.00094368642146534 + o[7]*(0.00060003561586052 + (-0.000095322787813974
         + o[1]*(8.8283690661692e-6 + 1.45389992595188e-15*o[9]))*tau1))) + o[
        5]*(-0.000283190801238040 + o[1]*(0.00060706301565874 + o[6]*(
        0.0189900682184190 + tau1*(0.032529748770505 + (0.0218417171754140 +
        0.000052838357969930*o[1])*tau1))));
      pi1 = 7.1000000000000 - pi;
      tau =  TSTAR1 /T;
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
  end DensityWater;

  model SaturationDensityDIPPR "saturation density using equation DIPPR 105"
  input SI.Temperature T;
  input Real[4] C;
  input SI.MolarMass MM;
  output SI.Density rho;

  Real Y;
  equation
  rho = C[1]/C[2]^Y * MM;
  Y= 1+ (1- T/C[3])^C[4];
  end SaturationDensityDIPPR;

  model SaturationDensityInfochem
    "saturation density using correlation from Infochem"
  input SI.Temperature T;
  input Real C[3];
  input SI.MolarMass MM;
  input SI.Temperature Tc;
  output SI.Density rho;
  Real tau;
  equation
  rho = (C[1] + C[2]*tau^C[3])*MM;
  tau = 1- T/Tc;
  end SaturationDensityInfochem;

  model TestActivity
   // SI.Temperature T=331.1;
   Real T=336.05;
   // Real x[4]={0.6, 0.233965-0.18, 0.178727, 0.162430};
   //  Real x[4]={0.424879, 0.233965, 0.178727, 0.162430};

  // Real x[:]={0.101, 0, 0.827, 0.072};
  // Real x[:]= { 0.34, 0, 0.209, 1-0.209-0.34};
   //Real x[:] = {0.762, 0, 0.136, 1-0.762-0.136};
  // Real x[:] = {0.616,0,0.205,1-0.616-0.205};
    Real x2[:] = {0.762,  0.136, 1-0.762-0.136};
  Real x[:] = {0.109,0, 0.442, 1-0.109-0.442};
  ThermalSeparation.Media.Methylacetatsynthese_Liq.ActivityCoefficient
      activityCoeff(                                                                 T=T,x_l=x);

  ThermalSeparation.Media.Correlations.ActivityCoefficient.Wilson wilsonHuss(nS=4, T=T, x_l=x, A={{0, -696.5, -31.19, 645.7}, {1123.1, 0, 2535.2, 237.5}, {813.2, -547.5, 0, 107.4}, {1917.2, 658, 469.55, 0}});

  ThermalSeparation.Media.Correlations.ActivityCoefficient.Wilson wilsonXu(nS=4, T=T, x_l=x,    A={{0, -696.5, -125.7, 887.4}, {1123.1, 0, 2.15, -65.8},   {891.1, 1.71, 0, 216.85}, {1966.7, 870.2, 468.6, 0}});

  ThermalSeparation.Media.Correlations.ActivityCoefficient.Wilson wilsonMartin(nS=3, T=T, x_l=x2,    A={{0, -45.55, 240.5}, {406.1, 0, -360.4},   {1013, 580.3,0}}, V_i={79.84, 44.44, 18.07});

  ThermalSeparation.Media.Correlations.ActivityCoefficient.NRTL nrtlXu(nS=4, T=T, x_l=x, alpha = {{0, 0.36, 1.029, 0.383},{0.36, 0, 0.3055, 0.2960},{1.0293, 0.3055,0, 0.2989}, {0.383, 0.36, 0.2989, 0}},
  g={{0,1218.87, 566.2, 879.1},{-635.9, 0, 0.577, -187.7},{456.9, -39.6, 0,-245.9},{1709, 981, 921,0}});

  ThermalSeparation.Media.Correlations.ActivityCoefficient.NRTL nrtlMartin(nS=3, T=T, x_l=x2, alpha = alphaMartin,
  g={{0, 168.9, 521.3},{177.6, 0, 333.2},{666.7, -121.7,0}});
  Real alphaMartin[3,3];
  parameter Real gMartin[3,3]={{0, 168.9, 521.3},{177.6, 0, 333.2},{666.7, -121.7,0}};
  equation
  for i in 1:3 loop
    for j in 1:3 loop

        alphaMartin[i,j]=if i==j then 0 else 1/(2+gMartin[i,j]*gMartin[j,i]);

      end for;
      end for;
  end TestActivity;

  model SatPressureForner
    parameter Integer n=4;
    input SI.Temperature T;
    output SI.Pressure p_sat[n];
  /*** saturation pressure ***/
    Real A_Ant[:] = {16.2685, 17.4068, 18.6073, 18.585};
    Real B_Ant[:] = {2665.5416, 3785.565, 3643.3113, 3984.9228};
    Real C_Ant[:] = {-53.424, -39.626, -33.4240, -39.724};

  equation
  for i in 1:n loop
    p_sat[i] = 133.3*exp(A_Ant[i] - B_Ant[i]/(T+C_Ant[i]));
  end for;

  end SatPressureForner;

  model SatPressureHuss
    parameter Integer n=4;
    input SI.Temperature T;
    output SI.Pressure p_sat[n];
  /*** saturation pressure ***/
    Real A_Ant[:] = {21.152, 22.1, 23.5, 23.22};
    Real B_Ant[:] = {-2662.78, -3654.6, -3643.3, -3835.2};
    Real C_Ant[:] = {-53.46, -45.4, -33.4, -45.3};

  equation
  for i in 1:n loop
    p_sat[i] = exp(A_Ant[i] + B_Ant[i]/(C_Ant[i]+T));
  end for;

  end SatPressureHuss;

  model SatPressureMartin
    parameter Integer n=3;
    input SI.Temperature T;
    output SI.Pressure p_sat[n];
  /*** saturation pressure ***/
    Real A_Ant[:] = {14.2, 15.8, 16.6};
    Real B_Ant[:] = {2665.5, 3242.8, 3984.9};
    Real C_Ant[:] = {-53.4, -49.55, -39.72};

  equation
  for i in 1:n loop
    p_sat[i] = exp(A_Ant[i] - B_Ant[i]/(C_Ant[i]+T))*1000;
  end for;

  end SatPressureMartin;

  model TestSatPressure
    SI.Temperature T=65+273;
  ThermalSeparation.Media.Methylacetatsynthese_Liq.SatPressureForner forner(T=T);
  ThermalSeparation.Media.Methylacetatsynthese_Liq.SatPressureHuss huss(T=T);
  ThermalSeparation.Media.Methylacetatsynthese_Liq.SatPressureMartin martin(T=T);

  end TestSatPressure;

redeclare model extends HenryCoefficient

end HenryCoefficient;

  model MolarSaturationDensityDIPPR
    "saturation density using equation DIPPR 105"
  input SI.Temperature T;
  input Real[4] C;
  input SI.MolarMass MM;
  output SI.Density rho;

  Real Y;
  equation
  rho = C[1]/C[2]^Y;
  Y= 1+ (1- T/C[3])^C[4];
  end MolarSaturationDensityDIPPR;

  model MolarSaturationDensityInfochem
    "saturation density using correlation from Infochem"
  input SI.Temperature T;
  input Real C[3];
  input SI.MolarMass MM;
  input SI.Temperature Tc;
  output SI.Density rho;
  Real tau;
  equation
  rho = (C[1] + C[2]*tau^C[3]);
  tau = 1- T/Tc;
  end MolarSaturationDensityInfochem;

  redeclare model extends DiffusionCoefficient
  /*** model which shall be used to calculate the diffusion coefficients in a binary mixture ***/
  model D_Molecules =
      ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.Molecule.Const_Molecule(D_const_liquid=1e-8);

  /*** model where the diffusion coefficients in a binary mixture are used to calculate the diffusion coefficients in a multicomponent mixture ***/
  ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.DiffMolecules
      diffCoeff(
      redeclare model D_Molecules = D_Molecules,
    T=T,
    p=p,
    x=x,
    nS=nSubstance,
    has_etaSubstance=has_etaSubstance,
    eta=eta);
  equation
  D = diffCoeff.D;
  end DiffusionCoefficient;

  redeclare replaceable model extends ThermodynamicFactor
     parameter Real r_UNI[:] = {2.8042, 2.2024, 1.4311, 0.92};
     parameter Real q_UNI[:] = {2.576, 2.072, 1.432, 1.4};
  //   parameter Real r_UNI[:] = {2.8042,1.9, 1.4311, 0.92};
  //   parameter Real q_UNI[:] = {2.576, 1.8, 1.432, 1.4};
    Real a_UNI[:,:] = {{1, 195.08,  620.11, 593.7},    {-76.614, 1,390.26, 422.38},        {  -86.02, 65.245, 1, -235.52},  {    -265.83, -98.12, -761.48, 1}};
    Real b_UNI[:,:] = {{1, 0.40148,-0.96715, 0.010143},{ -0.33388, 1, 0.97039, -0.051007}, { 0.13943, -2.0346, 1, 0.41274}, { 0.96295, -0.29355, 5.2418, 1}};
    Real c_UNI[:,:] = {{1, 0, 0, -2.1609e-3},          {  0, 1, -3.0613e-3, -2.4019e-4},   { 0, 3.157e-3, 1, -1.5415e-3},   {   2.0113e-4, -7.6741e-5, -4.3663e-3, 1}};
    Real u_UNI[Methylacetatsynthese_Liq.n,Methylacetatsynthese_Liq.n];

   ThermalSeparation.Media.Correlations.ThermodynamicFactorLiq.UNIQUAC
      thermodynamicFactorLiq(nS=n,T=T,x=x,q=q_UNI, r=r_UNI, lambda=u_UNI);
   Real omega = 0.5*tanh(0.01*(time-200))+0.5;
  equation
    /*** activity coefficient -ok ***/
  for i in 1:n loop
    for j in 1:n loop
      if i == j then
        u_UNI[i,j]=1;
      else

          u_UNI[i,j] = (a_UNI[i,j] + b_UNI[i,j]*T + c_UNI[i,j]*T^2)*Modelica.Constants.R;

      end if;
    end for;
        end for;
  Gamma=thermodynamicFactorLiq.Gamma;
  //Gamma=diagonal(ones(nSubstance));
  end ThermodynamicFactor;

         redeclare model extends ReactionRates
           Correlations.Reaction.ReactionRates.Langmuir langmuir(
             nS=4,
             nR=1,
             c=c,
             gamma=gamma,
             K_ads={4.15,3.15,5.64,5.24},
             x=propsLiq.x,
             T=propsLiq.T,
             MMX=MMX,
             rrc_forward(k_0=8.497e6*ones(nR), E=60.47e3*ones(nR)),
             rrc_reverse(k_0=6.127e5*ones(nR), E=63.73e3*ones(nR)),
             powerCoeffForward={{0,1,1,0}},
             powerCoeffReverse={{1,0,0,1}});

           input ThermodynamicProperties propsLiq;
           //K_ads aus Forner
         equation
           r = langmuir.r;
           k_forward= langmuir.k_forward;
         end ReactionRates;

         redeclare model extends EquilibriumConstant

               /*** reaction rate coefficient ***/
           replaceable model RRC_forward =
        ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.Arrhenius
                                                                      constrainedby
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.BaseReactionRateCoeff
                                                                                                          annotation(choicesAllMatching=true);

           RRC_forward rrc_forward(nR=nR, T= propsLiq.T, k_0=8.497e6*ones(nR), E=60.47e3*ones(nR));

          replaceable model RRC_reverse =
        ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.Arrhenius
                                                                      constrainedby
      ThermalSeparation.Media.Correlations.Reaction.ReactionRates.ReactionRateCoeff.BaseReactionRateCoeff
                                                                                                          annotation(choicesAllMatching=true);

           RRC_reverse rrc_reverse(nR=nR, T = propsLiq.T, k_0=6.127e5*ones(nR), E=63.73e3*ones(nR));

             Real k_forward[nR]= rrc_forward.k "reaction rate coefficient forward reaction";
           Real k_reverse[nR]=rrc_reverse.k "reaction rate coefficient reverse reaction";

         Correlations.Reaction.EquilibriumConstant.Langmuir eqConst(
             nS=nSubstance,
             nR=1,
             c=c,
             gamma=gamma,
                 K_ads={4.15,3.15,5.64,5.24},
             x=propsLiq.x,
             T=propsLiq.T,
             MMX=MMX,
             powerCoeffForward={{0,1,1,0}},
             powerCoeffReverse={{1,0,0,1}});

               ThermodynamicProperties propsLiq;

         equation
           K_eq = k_forward./k_reverse;
           K_eq_MWG = eqConst.K_eq_MWG;
         end EquilibriumConstant;
end Methylacetatsynthese_Liq;
