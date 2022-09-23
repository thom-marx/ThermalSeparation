within ThermalSeparation.Media.IdealGasMixtures;
package H2O_CO2 "H2O, CO2 for Siemens"

constant Real f_psat = 0.907;
constant Boolean psat_Antoine = false;
     constant Integer henry[:] = {2};
     constant SI.Pressure henry_H[:] = {0,1.63e9} "Henry coefficient";
     constant Real henry_C[:]={0,2400}
    "constant to calculate temperature dependency";
     constant SI.Temperature henry_T_norm=298 "norm temperature";

  import
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices.ReferenceEnthalpy;

  extends ThermalSeparation.Media.IdealGasMixtures.BaseClasses.BaseIdealGasMixture(
                                                                             reference_X={0.92,0.08},
                              nSubstance=2,  V = {13.1,  26.7}, eq_Tsonopoulos = {6, 1},  sigma = {2.641,  3.941}, epsilon_k = {809.1, 195.2}, nS=2,
    SpecificEnthalpy(start=if referenceChoice == ReferenceEnthalpy.ZeroAt0K then
                3e5 else if referenceChoice == ReferenceEnthalpy.UserDefined then
                h_offset else 0, nominal=1.0e5),
    Density(start=2, nominal=10),
    AbsolutePressure(start= 10e5, nominal=10e5),
    Temperature(start=500, nominal=500),
    fluidConstants={ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.H2O,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.CO2},
                               data=
      {ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.H2O,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.CO2},
                                  R_const=data.R_s);

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

  protected
   parameter Units.ConversionDebyeCm conversionDebye = 3.33564e-30
      "to convert from debye to C*m";
  equation
  for i in 1:nX loop
   X[i] = max(1e-5,min(1,x[i]*MMX[i]/MM));
  end for;

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

   output SI.Temperature T_sat_water;
  protected
   Real hX[nSubstance];
    /***Berechnung der Sttigungstemperatur von Wasser beim Partialdruck des Wasserdampfes***/
    Real pi "dimensionless pressure";
    Real[20] o "vector of auxiliary variables";

    /*** evaporation enthalpy of water ***/
    SI.SpecificEnthalpy r_water=-2462.5*T + 3177.8e3;

    /*** heat capacity ***/
      parameter SI.SpecificHeatCapacity cp_water = 4200;

  equation
        /***calculation of the saturation temperature of water at the partial pressure of the water vapour***/
      //Vorkehrungen treffen, falls der Wasseranteil null ist (x[1] = 0)
      pi =  max(1e-4,p*x[1])*1e-6/f_psat;
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

          hX[1] = (cp_water*(T_sat_water-T0)  + 1875*(T-T_sat_water)) *MMX[1];

      hX[2] = 864*MMX[2] * (T-T0);

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;

  redeclare replaceable model extends EvaporationEnthalpy

  protected
    SI.SpecificEnthalpy r_water=-2462.5*T + 3177.8e3;
  equation
      /*** enthalpy ***/
    h[1] = MMX[1]*r_water;
    h[2] = 0;

  annotation (Icon(graphics={
                        Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end EvaporationEnthalpy;
end H2O_CO2;
