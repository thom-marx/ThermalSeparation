within ThermalSeparation.Media;
package C2H5OH_Water_Liq "2 components: C2H5OH, H2O"
  extends BaseMediumLiquid(has_etaSubstance={true,true}, nSubstance=2, Tcrit= {513.9, 647.3}, pcrit= {61.48e5, 2.2064e7}, Vcrit = {167.00e-6, 55.95e-6},omega = { 0.037,
 0.022}, mu = {1.7, 1.8}, MMX= {0.046,0.018}, eq_Tsonopoulos = {4, 6}, henry=fill(false,nSubstance));

  constant Integer n = nSubstance "number of components";

    constant SI.MolarMass MM_substances[n] = MMX;

  redeclare replaceable model extends BaseProperties
    constant SI.SpecificHeatCapacity RH2O=461.526
      "specific gas constant of water vapour";
    constant SI.MolarMass MH2O=0.01801528 "molar weight of water";
    constant SI.Temperature TSTAR1=1386.0
      "normalization temperature for region 1 IF97";
    constant SI.Pressure PSTAR1=16.53e6
      "normalization pressure for region 1 IF97";
    constant SI.Temperature TSTAR2=540.0
      "normalization temperature for region 2 IF97";
    constant SI.Pressure PSTAR2=1.0e6
      "normalization pressure for region 2 IF97";
    constant SI.Temperature TSTAR5=1000.0
      "normalization temperature for region 5 IF97";
    constant SI.Pressure PSTAR5=1.0e6
      "normalization pressure for region 5 IF97";
    constant SI.SpecificEnthalpy HSTAR1=2.5e6
      "normalization specific enthalpy for region 1 IF97";
    constant Real IPSTAR=1.0e-6
      "normalization pressure for inverse function in region 2 IF97";

  constant SI.DynamicViscosity sigma_EtOH=0.02255;// bei 20 °C N/m
  constant SI.DynamicViscosity sigma_H2O=0.0728;// bei 20 °C N/m

  Real MM_min;
  Real MM_var;
  SI.MassFraction w[n](start=fill(0.03,n));

     CalcSpecificEnthalpy calcSpecificEnthalpy(T0=T0, p=p, T=T, x=x);
  //Real lambda;
  Real lambda_EtOH = 0.155;
  Real lambda_H2O= 0.66;
  /***density and enthalpy ***/
  DensityWater densityWater(p=1e5,T=T,scale=1);
  Real d_molar_water=d2/MMX[2];
  Real d_molar_ethanol=d1/MMX[1];

  protected
    Real d2 = densityWater.d;
    Real d1;

  equation
    eta_comp= fill(eta,n);
    v = x[1]/d_molar_ethanol + x[2]/d_molar_water;
        for j in 1:n loop
    w[j] =MMX[j]/MM*x[j];
    end for;
    eta=0.001002;
    sigma=min(0.0728,max(0.02255,sigma_EtOH*w[1]+sigma_H2O*w[2]));

   d=1/((1/max(1e-5,d1))*max(1e-5,w[1])+(1/max(1e-5,d2))*max(1e-5,w[2]));
    MM_min=min(MMX);
    MM= MM_var;
    MM_var = max(MM_min,sum(x[j]*MMX[j] for j in 1:n));

    h = calcSpecificEnthalpy.h;//h_Jkg * MM;
    u = h - p/d*MM;
    p_sat[1]=10^(8.11220-1592.864/((T-273.15)+226.184))*133.322; // Saturationpressure for Ethanol Perry's Chemical Engineers 13-14 in Pa (T=20-93°C)
    p_sat[2]=10^(8.07131-1730.63/((T-273.15)+233.426))*133.322; // Saturationpressure for Water Perry's Chemical Engineers 13-14 in Pa (T=1-100°C)

            /*** density ***/
            d1=761.2; // bei 50 °C aus VDI Wärmeatlas
  //     d2 = p/( RH2O * T * pi * gpi);
  //     pi = p/ PSTAR1;
  //     gpi = pi1*(pi1*(o[10]*(0.000095038934535162 + o[2]*(
  //       8.4812393955936e-6 + 2.55615384360309e-9*o[6])) + pi1*(o[12]*(
  //       8.9701127632000e-6 + (2.60684891582404e-6 + 5.7366919751696e-13*o[13])
  //       *o[7]) + pi1*(2.02584984300585e-6*o[14] + o[16]*((1.01874413933128e-8
  //        + 1.39398969845072e-9*o[11])*o[36] + o[19]*(1.44400475720615e-17*o[
  //       34] + o[15]*(-3.3300108005598e-19*o[32] + o[20]*(-7.6373766822106e-22
  //       *o[30] + pi1*(3.5842867920213e-22*o[28] + pi1*(-5.6507093202352e-23*o[
  //       26] + 2.99318679335866e-24*o[24]*pi1))))))))) + o[8]*(
  //       0.00094368642146534 + o[7]*(0.00060003561586052 + (-0.000095322787813974
  //        + o[1]*(8.8283690661692e-6 + 1.45389992595188e-15*o[9]))*tau1))) + o[
  //       5]*(-0.000283190801238040 + o[1]*(0.00060706301565874 + o[6]*(
  //       0.0189900682184190 + tau1*(0.032529748770505 + (0.0218417171754140 +
  //       0.000052838357969930*o[1])*tau1))));
  //     pi1 = 7.1000000000000 - pi;
  //     tau =  TSTAR1 /T;
  //     tau1 = -1.22200000000000 + tau;
  //       o[1] =  tau1*tau1;
  //   o[2] = o[1]*o[1];
  //   o[3] = o[2]*o[2];
  //   o[4] = o[3]*tau1;
  //   o[5] = 1/o[4];
  //   o[6] = o[1]*o[2];
  //   o[7] = o[1]*tau1;
  //   o[8] = 1/o[7];
  //   o[9] = o[1]*o[2]*o[3];
  //   o[10] = 1/o[2];
  //   o[11] = o[2]*tau1;
  //   o[12] = 1/o[11];
  //   o[13] = o[2]*o[3];
  //   o[14] = 1/o[3];
  //   o[15] = pi1*pi1;
  //   o[16] = o[15]*pi1;
  //   o[17] = o[15]*o[15];
  //   o[18] = o[17]*o[17];
  //   o[19] = o[17]*o[18]*pi1;
  //   o[20] = o[15]*o[17];
  //   o[21] = o[3]*o[3];
  //   o[22] = o[21]*o[21];
  //   o[23] = o[22]*o[3]*tau1;
  //   o[24] = 1/o[23];
  //   o[25] = o[22]*o[3];
  //   o[26] = 1/o[25];
  //   o[27] = o[1]*o[2]*o[22]*tau1;
  //   o[28] = 1/o[27];
  //   o[29] = o[1]*o[2]*o[22];
  //   o[30] = 1/o[29];
  //   o[31] = o[1]*o[2]*o[21]*o[3]*tau1;
  //   o[32] = 1/o[31];
  //   o[33] = o[2]*o[21]*o[3]*tau1;
  //   o[34] = 1/o[33];
  //   o[35] = o[1]*o[3]*tau1;
  //   o[36] = 1/o[35];
  //   o[37] = o[1]*o[3];
  //   o[38] = 1/o[37];
  //   o[39] = 1/o[6];
  //   o[40] = o[1]*o[22]*o[3];
  //   o[41] = 1/o[40];
  //   o[42] = 1/o[22];
  //   o[43] = o[1]*o[2]*o[21]*o[3];
  //   o[44] = 1/o[43];
  //   o[45] = 1/o[13];

      lambda=x[1]*lambda_EtOH+x[2]*lambda_H2O;
      cp=4000;

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end BaseProperties;

redeclare replaceable model extends ActivityCoefficient

  replaceable model ActivityCoeff =
      ThermalSeparation.Media.Correlations.ActivityCoefficient.Ideal                                     constrainedby
      ThermalSeparation.Media.Correlations.ActivityCoefficient.BaseActivityCoefficient
    annotation (choicesAllMatching=true, Dialog(tab="General", group=
          "Thermodynamic Equilibrium"));
   ActivityCoeff activityCoeff(nS=nSubstance, T=T, x_l = x_l);

equation
  gamma = activityCoeff.gamma;

end ActivityCoefficient;

redeclare model extends FugacityCoefficient
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

end FugacityCoefficient;

  model Testmodell
    parameter Integer n=1;
    parameter Integer nS=2;
  parameter SI.Temperature T = 373.15;
  parameter SI.Pressure p = 100000;
  parameter SI.MoleFraction x[nS] = {0.3,0.7};//,0.05}};

  BaseProperties BP(
      T=T,
      p=p,
      x=x);
  end Testmodell;

  redeclare model extends CalcSpecificEnthalpy
  /*** enthalpy ***/

  constant SI.MolarHeatCapacity cp_l_EtOH = 2443.48 *MMX[1];
  constant SI.MolarHeatCapacity cp_l_H2O = 4187*MMX[2];
   ThermalSeparation.Units.MolarEnthalpy h_l_EtOH;
   ThermalSeparation.Units.MolarEnthalpy h_l_H2O;

    SI.MolarMass MM = max(1e-7,sum(x.*MM_substances));

  equation
        /*** enthalpy ***/

      h_l_EtOH=cp_l_EtOH*(T-T0);
    h_l_H2O=cp_l_H2O*(T-T0);
      h = x[1]*h_l_EtOH+x[2]*h_l_H2O;

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;

redeclare model extends HenryCoefficient

end HenryCoefficient;

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

  redeclare model extends DiffusionCoefficient(diluted=true)
    /*** model which shall be used to calculate the diffusion coefficients in a binary mixture ***/
    model D_Molecules =
        ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.Molecule.Const_Molecule;

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
    D_diluted[1]=diffCoeff.D0[1,2];
    D_diluted[2]=diffCoeff.D0[2,1];
  end DiffusionCoefficient;

  redeclare replaceable model extends ThermodynamicFactor
  equation
    Gamma= diagonal(ones(nSubstance));
  end ThermodynamicFactor;
end C2H5OH_Water_Liq;
