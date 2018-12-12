within ThermalSeparation.Media.WaterBasedLiquid.BaseClasses;
partial package PartialWaterBasedReaction

  extends
    ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.PartialMediumReaction;
  constant Boolean smoothModel
    "true if the (derived) model should not generate state events";
  constant Boolean onePhase
    "true if the (derived) model should never be called with two-phase inputs";

  record FluidLimits "validity limits for fluid model"
    import ThermalSeparation;
    extends Modelica.Icons.Record;
    ThermalSeparation.Media.Types.Temperature TMIN "minimum temperature";
    ThermalSeparation.Media.Types.Temperature TMAX "maximum temperature";
    ThermalSeparation.Media.Types.Density DMIN "minimum density";
    ThermalSeparation.Media.Types.Density DMAX "maximum density";
    ThermalSeparation.Media.Types.AbsolutePressure PMIN "minimum pressure";
    ThermalSeparation.Media.Types.AbsolutePressure PMAX "maximum pressure";
    ThermalSeparation.Media.Types.SpecificEnthalpy HMIN "minimum enthalpy";
    ThermalSeparation.Media.Types.SpecificEnthalpy HMAX "maximum enthalpy";
    ThermalSeparation.Media.Types.SpecificEntropy SMIN "minimum entropy";
    ThermalSeparation.Media.Types.SpecificEntropy SMAX "maximum entropy";
    annotation(Documentation(
        info="<html>
          <p>The minimum pressure mostly applies to the liquid state only.
          The minimum density is also arbitrary, but is reasonable for techical
          applications to limit iterations in non-linear systems. The limits in
          enthalpy and entropy are used as safeguards in inverse iterations.</p>
          </html>"));
  end FluidLimits;

  redeclare replaceable model extends BaseProperties
    "Base properties (p, d, T, h, u, R, MM, sat) of two phase medium"
    SaturationProperties sat "Saturation properties at the medium pressure";

    ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.data
      data;
    DensityWater densityWater(p=p,T=T);

  protected
    SI.Density d_water;
  equation
    eta_comp= fill(eta,nSubstance);
    d_water = densityWater.d;
    v = MM/d;
      lambda = thermalConductivity(state);
      cp=specificHeatCapacityCp(state);
    annotation(Documentation(info="<html></html>"));
  end BaseProperties;

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

    end ActivityCoefficient;

        redeclare replaceable model extends FugacityCoefficient "calculates the fugacity coefficient at the saturation point for each component in the liquid phase"

            replaceable model FugacityCoeff =
        ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient.SaturationFugacitycoefficient
                                                                                                          constrainedby
      ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient.BaseFugacityCoefficient
            annotation (choicesAllMatching=true, Dialog(tab="General", group=
                  "Thermodynamic Equilibrium"));

             FugacityCoeff fugacityCoeff(nS=nSubstance, T=T,p=p, p_sat=p_sat, Tcrit=Tcrit, pcrit=pcrit, omega=omega, NoOfEq=eq_Tsonopoulos, mu=mu);

        equation
          phi_sat = fugacityCoeff.phi_sat;
          annotation (Icon(graphics={
                            Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}));
        end FugacityCoefficient;

  model DensityWater
    input SI.Pressure p;
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

      d =p/( RH2O * T * pi * gpi);
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

  redeclare replaceable model extends ThermodynamicFactor
  equation
  Gamma=diagonal(ones(nSubstance));
  end ThermodynamicFactor;
end PartialWaterBasedReaction;
