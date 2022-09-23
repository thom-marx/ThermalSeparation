within ThermalSeparation.Media;
package H2O_CO2_Vap "H2O,CO2: CO2 separation with MEA"
     constant Integer henry[:] = {2};
                                     //1.63e9
     constant Modelica.Units.SI.Pressure henry_H[:]={0,1.63e9}
  "Henry coefficient";
     constant Real henry_C[:]={0,2400}
  "constant to calculate temperature dependency";
     constant Modelica.Units.SI.Temperature henry_T_norm=298 "norm temperature";

import
  ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices.ReferenceEnthalpy;


  extends ThermalSeparation.Media.IdealGasMixtures.BaseClasses.BaseIdealGasMixture(
    reference_X={0.92,0.08},
    nSubstance=2,
    delta_hv_medium=true,
    V={13.1,26.7},
    eq_Tsonopoulos={6,1},
    sigma={2.641,3.941},
    epsilon_k={809.1,195.2},
    nS=2,
    SpecificEnthalpy(start=if referenceChoice == ReferenceEnthalpy.ZeroAt0K then
                3e5 else if referenceChoice == ReferenceEnthalpy.UserDefined then
                h_offset else 0, nominal=1.0e5),
    Density(start=2, nominal=10),
    AbsolutePressure(start=10e5, nominal=10e5),
    Temperature(start=500, nominal=500),
    fluidConstants={ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.H2O,ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.CO2},
    data={ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.H2O,ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.CO2},
    R_const=data.R);


  redeclare replaceable model extends BaseProperties(
    T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    Xi(each stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default))
  import
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasNasa;
    CalcSpecificEnthalpy calcSpecificEnthalpy(
      T0=T0,
      T=T,
      x=x,
      p=p);

  Real p_start(start=1.0e5);

    //   /*** Enthalpy ***/
    //   Modelica.Units.SI.SpecificHeatCapacity cp;
    //   // parameter SI.Temperature T0 = 293.15 "reference temperature";
    //
    //   /*** thermal conductivity ***/
    //     Modelica.Units.SI.ThermalConductivity lambda=0.02;
protected
   parameter ThermalSeparation.Units.ConversionDebyeCm conversionDebye=
      3.33564e-30 "to convert from debye to C*m";
  equation
  for i in 1:nX loop
   X[i] = max(1e-5,min(1,x[i]*MMX[i]/MM));
  end for;
  //cp=specificHeatCapacityCp(state);

  p_start=state.p;
    assert(T >= 200 and T <= 6000, "
Temperature T (="   + String(T) + " K = 200 K) is not in the allowed range
200 K <= T <= 6000 K
required from medium model \""   + mediumName + "\".");

    //bei der Gleichung MM=molarMass(state): Index-Problem!
  // MM = max(min(MMX),sum(x.*MMX));//molarMass(state);
  d=sum(c.*MMX);

    h =calcSpecificEnthalpy.h;

    R = data.R_s*X;

    u = h-R*T*MM;
    d = p/(R*T);

    annotation (structurallyIncomplete);
  end BaseProperties;

  redeclare replaceable model extends CalcSpecificEnthalpy
  constant Real f_psat = 1;
  constant Boolean psat_Antoine = false;
   output Modelica.Units.SI.Temperature T_sat_water;
protected
   Modelica.Units.SI.SpecificEnthalpy hX[nSubstance];
    /***Berechnung der Sättigungstemperatur von Wasser beim Partialdruck des Wasserdampfes***/
    Real pi "dimensionless pressure";
    Real[20] o "vector of auxiliary variables";

    /*** evaporation enthalpy of water ***/
    Modelica.Units.SI.SpecificEnthalpy r_water;
    Real tau;
    Real chvap[5];

    /*** heat capacity ***/
      parameter Modelica.Units.SI.SpecificHeatCapacity cp_water=4200;
      Real cp_H2O_vap;
      Real cp_CO2_vap;
      Real cp_coeff_H2O[5];
      Real cp_coeff_CO2[5];

  equation
        /***Berechnung der Sättigungstemperatur von Wasser beim Partialdruck des Wasserdampfes***/
      //Vorkehrungen treffen, falls der Wasseranteil null ist (x[1] = 0)
      pi =  max(1e-4,p*x[1])*1e-6/f_psat;//min(p,data.PCRIT)*data.IPSTAR;
    o[1] =  pi^0.25;
    o[2] =  -3.2325550322333e6*o[1];
    o[3] =  max(1e-10,pi)^0.5;
    o[4] =  -724213.16703206*o[3];
    o[5] =  405113.40542057 + o[2] + o[4];
    o[6] =  -17.0738469400920*o[1];
    o[7] =  14.9151086135300 + o[3] + o[6];
    o[8] =  -4.0*o[5]*o[7];
    o[9] =  12020.8247024700*o[1];
    o[10] =  1167.05214527670*o[3];
    o[11] =  -4823.2657361591 + o[10] + o[9];
    o[12] =  o[11]*o[11];
    o[13] =  o[12] + o[8];
    o[14] =  max(1e-8,o[13])^0.5;
    o[15] =  -o[14];
    o[16] =  -12020.8247024700*o[1];
    o[17] =  -1167.05214527670*o[3];
    o[18] =  4823.2657361591 + o[15] + o[16] + o[17];
    o[19] =  1/o[18];
    o[20] =  2.0*o[19]*o[5];

    if psat_Antoine then
        /*** Antoine-Gleichung von Siemens ***/
    T_sat_water = 5342.9/(24.9022 - Modelica.Math.log(max(1e-7,p*x[1]))-21.2826);
    else
      T_sat_water = max(200,0.5*(650.17534844798 + o[20] - (-4.0*(-0.238555575678490 +
       1300.35069689596*o[19]*o[5]) + (650.17534844798 + o[20])^2.0)^0.5));
  end if;

        h = x[1]*hX[1] +  x[2]*hX[2];

          //hX[1] = (cp_water*(T_sat_water-T0)  + 1875*(T-T_sat_water)) *MMX[1] + MMX[1] * r_water;
          //hX[1] = (cp_water*(T-T0)  + 1875*(T-T)) *MMX[1] + MMX[1] * r_water; // hier T_sat_water = T da Gleichgewicht
          hX[1] = cp_H2O_vap*(T-T0) + r_water;  // Oexmann doesn't consider liquid part of gas enthalpy

      hX[2] = cp_CO2_vap * (T-T0);

        tau = 1 - (T/647.3);
        r_water = 8.31451 * 647.3 *
                        (chvap[1] * max(1e-7,tau)^(1.0/3.0) +
                                        chvap[2] * max(1e-7,tau)^(2.0/3.0)  +
                                        chvap[3] * max(1e-7,tau) +
                                        chvap[4] * max(1e-7,tau)^(2.0)  +
                                        chvap[5] * max(1e-7,tau)^(6.0));

        // coefficients for water vaporisation
        chvap ={5.6297, 13.962, -11.673, 2.1784, -0.31666};

        cp_H2O_vap = cp_coeff_H2O[1] + cp_coeff_H2O[2] *((cp_coeff_H2O[3]/T)/Modelica.Math.sinh(cp_coeff_H2O[3]/T))^2 + cp_coeff_H2O[4] * ((cp_coeff_H2O[5]/T)/Modelica.Math.cosh(cp_coeff_H2O[5]/T))^2;
        cp_CO2_vap = cp_coeff_CO2[1] + cp_coeff_CO2[2] *((cp_coeff_CO2[3]/T)/Modelica.Math.sinh(cp_coeff_CO2[3]/T))^2 + cp_coeff_CO2[4] * ((cp_coeff_CO2[5]/T)/Modelica.Math.cosh(cp_coeff_CO2[5]/T))^2;

        cp_coeff_H2O ={33.363,26.79,2610.5,8.896,1169};
        cp_coeff_CO2 ={29.37,34.54,1428,26.4,588};
    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;

  redeclare replaceable model extends EvaporationEnthalpy

protected
    Modelica.Units.SI.SpecificEnthalpy r_water;
    Real tau;
    Real chvap[5];
  equation
      /*** enthalpy ***/
    h[1] = r_water;
    h[2] = 0;

        tau = 1 - (T/647.3);
        r_water = 8.31451 * 647.3 *
                        (chvap[1] * tau^(1.0/3.0) +
                                        chvap[2] * tau^(2.0/3.0)  +
                                        chvap[3] * tau +
                                        chvap[4] * tau^(2.0)  +
                                        chvap[5] * tau^(6.0));

        // coefficients for water vaporisation
        chvap ={5.6297, 13.962, -11.673, 2.1784, -0.31666};

  annotation (Icon(graphics={
                        Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end EvaporationEnthalpy;
end H2O_CO2_Vap;
