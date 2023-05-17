within ThermalSeparation.Media.IdealGasMixtures;
package N2_H2O_O2_CO2 "4 components: N2, H2O, O2, CO2"

       constant Integer henry[:] = {1,3,4};
                 constant SI.Pressure henry_H[:] = {8.51e9,0,4.25e9,1.63e8}
    "Henry coefficient";
     constant Real henry_C[:]={1300,0,1500,2400}
    "constant to calculate temperature dependency";
     constant SI.Temperature henry_T_norm=298 "norm temperature";

  import
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices.ReferenceEnthalpy;
  extends ThermalSeparation.Media.IdealGasMixtures.BaseClasses.BaseIdealGasMixture(
                                                                             delta_hv_medium = true,
                              nS=4,nSubstance=4,V = {18.5, 13.1, 16.3, 26.7}, eq_Tsonopoulos = {1, 6, 1, 1}, sigma = {3.798, 2.641, 3.467, 3.941},  epsilon_k={71.4, 809.1,106.7, 195.2},
    SpecificEnthalpy(start=if referenceChoice == ReferenceEnthalpy.ZeroAt0K then
                3e5 else if referenceChoice == ReferenceEnthalpy.UserDefined then
                h_offset else 0, nominal=1.0e5),
    Density(start=2, nominal=10),
    AbsolutePressure(start= 10e5, nominal=10e5),
    Temperature(start=500, nominal=500),
    fluidConstants={ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.N2,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.H2O,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.O2,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.CO2},
                                                                          data=
      {ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.N2,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.H2O,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.O2,
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.CO2},R_const=data.R_s);

 redeclare replaceable model extends BaseProperties(
    T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    Xi(each stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default))
    import
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasNasa;
          CalcSpecificEnthalpy calcSpecificEnthalpy(p=p,x=x,T=T);

  Real p_start(start=1.6e5);

  /***Berechnung der Sättigungstemperatur von Wasser beim Partialdruck des Wasserdampfes***/
  protected
    Real pi "dimensionless pressure";
    Real[20] o "vector of auxiliary variables";
    SI.Temperature T_sat_water;

  protected
   parameter Real conversionDebye = 3.33564e-30 "to convert from debye to C*m";
   SI.SpecificEnthalpy h_Jkg;
   SI.SpecificEnthalpy u_Jkg;
 equation
  for i in 1:nX loop
   X[i] = max(1e-5,min(1,x[i]*MMX[i]/MM));
  end for;
   ///// Tcrit = test();
  p_start=state.p;
    assert(T >= 200 and T <= 6000, "
Temperature T (="   + String(T) + " K = 200 K) is not in the allowed range
200 K <= T <= 6000 K
required from medium model \""   + mediumName + "\".");

      /***Berechnung der Sättigungstemperatur von Wasser beim Partialdruck des Wasserdampfes***/
      //Vorkehrungen treffen, falls der Wasseranteil null ist (x[1] = 0)
      pi =  max(1e-4,p*x[2])*1e-6;//min(p,data.PCRIT)*data.IPSTAR;
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
    T_sat_water =  max(300,0.5*(650.17534844798 + o[20] - (-4.0*(-0.238555575678490 +
      1300.35069689596*o[19]*o[5]) + (650.17534844798 + o[20])^2.0)^0.5));

    //bei der Gleichung MM=molarMass(state): Index-Problem!
   MM = max(1e-6,sum(x.*MMX));//molarMass(state);
    h_Jkg = X[2]*(4120*(T_sat_water-0) + 2382e3 + 1875*(T-T_sat_water))+ (X[3]*925 + X[4]*864  + X[1] * 1041)*(T-0);//h_TX(T, X);
    h = calcSpecificEnthalpy.h;
    R = data.R_s*X;
    u_Jkg = h_Jkg - R*T;
    u = u_Jkg * MM;
    d = p/(R*T);

    annotation (structurallyIncomplete);
 end BaseProperties;

  redeclare replaceable model extends CalcSpecificEnthalpy

  Real hX[nX];

  protected
    Real pi "dimensionless pressure";
    Real[20] o "vector of auxiliary variables";
    SI.Temperature T_sat_water;

  equation
        /***calculation of the saturation temperature of water at the partial pressure of the water vapour***/
          //Vorkehrungen treffen, falls der Wasseranteil null ist (x[1] = 0)
      pi =  max(1e-4,p*x[2])*1e-6;
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
    T_sat_water =  max(300,0.5*(650.17534844798 + o[20] - (-4.0*(-0.238555575678490 +
      1300.35069689596*o[19]*o[5]) + (650.17534844798 + o[20])^2.0)^0.5));

      h = sum(x.*hX);

      hX[2] = (4120*(T_sat_water-0) + 2382e3 + 1875*(T-T_sat_water))*MMX[2];
      hX[3] = 925*MMX[3] * (T-0);
      hX[4] = 864*MMX[4] * (T- 0);
      hX[1] = 1041*MMX[1]* (T-0);

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
      h[1]=0;
    h[2] = MMX[2]*r_water;
  for i in 3:nX loop
    h[i] = 0;
  end for;

  end EvaporationEnthalpy;
end N2_H2O_O2_CO2;
