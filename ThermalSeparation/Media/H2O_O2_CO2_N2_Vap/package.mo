within ThermalSeparation.Media;
package H2O_O2_CO2_N2_Vap "H2O,O2,CO2,N2: CO2 separation with MEA"
  constant Real f_psat = 0.907;
constant Boolean psat_Antoine = false;
     constant Integer henry[:] = {3};
                                     //1.63e9
     constant Modelica.SIunits.Pressure henry_H[:]={0,0,1.63e8,0}
  "Henry coefficient";
     constant Real henry_C[:]={0,1500,2400,1300}
  "constant to calculate temperature dependency";
     constant Modelica.SIunits.Temperature henry_T_norm=298 "norm temperature";

  import
  ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices.ReferenceEnthalpy;


  extends
  ThermalSeparation.Media.IdealGasMixtures.BaseClasses.BaseIdealGasMixture(
                              nS=4,nSubstance=4, V = {13.1, 16.3, 26.7, 18.5}, eq_Tsonopoulos = {6, 1, 1, 1},  sigma = {2.641, 3.467, 3.941, 3.798}, epsilon_k = {809.1,106.7, 195.2, 71.4},
    SpecificEnthalpy(start=if referenceChoice == ReferenceEnthalpy.ZeroAt0K then
                3e5 else if referenceChoice == ReferenceEnthalpy.UserDefined then
                h_offset else 0, nominal=1.0e5),
    Density(start=2, nominal=10),
    AbsolutePressure(start= 10e5, nominal=10e5),
    Temperature(start=500, nominal=500),
    fluidConstants={ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.H2O,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.O2,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.CO2,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.N2},             data=
      {ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.H2O,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.O2,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.CO2,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.N2},                            R_const=data.R);


 redeclare replaceable model extends BaseProperties(
    T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    Xi(each stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default))
    import
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasNasa;

    CalcSpecificEnthalpy calcSpecificEnthalpy(T0=T0,T=T, x=x, p=p);

  Real p_start(start=1.6e5);

  //   /*** Enthalpy ***/
  //   Modelica.SIunits.SpecificHeatCapacity cp;
  //
  //   /*** thermal conductivity ***/
  //     Modelica.SIunits.ThermalConductivity lambda=0.02;
protected
   parameter Real conversionDebye = 3.33564e-30 "to convert from debye to C*m";
 equation
  for i in 1:nX loop
   X[i] = max(1e-5,min(1,x[i]*MMX[i]/MM));
  end for;
 // cp=specificHeatCapacityCp(state);

  p_start=state.p;
    assert(T >= 200 and T <= 6000, "
Temperature T (="   + String(T) + " K = 200 K) is not in the allowed range
200 K <= T <= 6000 K
required from medium model \""   + mediumName + "\".");

    //bei der Gleichung MM=molarMass(state): Index-Problem!
  // MM = max(1e-6,sum(x.*MMX));//molarMass(state);
  d=sum(c.*MMX);

    h =calcSpecificEnthalpy.h;

    R = data.R*X;

    u = h-R*T*MM;

    d =   p/(R*T);

    annotation (structurallyIncomplete);
 end BaseProperties;


  redeclare replaceable model extends CalcSpecificEnthalpy

   output Modelica.SIunits.Temperature T_sat_water;

    /***Berechnung der Sättigungstemperatur von Wasser beim Partialdruck des Wasserdampfes***/

protected
   Modelica.SIunits.SpecificEnthalpy hX[4];
      Real pi "dimensionless pressure";
    Real[20] o "vector of auxiliary variables";

    /*** evaporation enthalpy of water ***/
    Modelica.SIunits.SpecificEnthalpy r_water=-2462.5*T + 3177.8e3;

    /*** heat capacity ***/
      parameter Modelica.SIunits.SpecificHeatCapacity cp_water=4200;

  equation
        /***Berechnung der Sättigungstemperatur von Wasser beim Partialdruck des Wasserdampfes***/
      //Vorkehrungen treffen, falls der Wasseranteil null ist (x[1] = 0)
    /*** Antoine-Gleichung von Siemens ***/

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
            T_sat_water = 5342.9/(24.9022 - Modelica.Math.log(max(1e-7,p*x[1]))-21.2826);
            else
      T_sat_water =  max(200,0.5*(650.17534844798 + o[20] - (-4.0*(-0.238555575678490 +
        1300.35069689596*o[19]*o[5]) + (650.17534844798 + o[20])^2.0)^0.5));
   end if;

        h = x[1]*hX[1] + x[3]*hX[3] + x[2]*hX[2] + x[4] * hX[4];

          hX[1] = (1875*(T-T0)) *MMX[1];
      hX[3] = 864*MMX[3] * (T-T0);
      hX[2] = 925*MMX[2] * (T- T0);
      hX[4] = 1041*MMX[4] * (T- T0);

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;


  redeclare replaceable model extends EvaporationEnthalpy

protected
    Modelica.SIunits.SpecificEnthalpy r_water=-2462.5*T + 3177.8e3;
  equation
      /*** enthalpy ***/
    h[1] = MMX[1]*r_water;
    h[2] = 0;
    h[3] = 0;
    h[4] = 0;

  annotation (Icon(graphics={
                        Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end EvaporationEnthalpy;
end H2O_O2_CO2_N2_Vap;
