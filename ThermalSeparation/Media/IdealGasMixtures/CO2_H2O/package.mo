within ThermalSeparation.Media.IdealGasMixtures;
package CO2_H2O "CO2, H2O"

     constant Integer henry[:] = {4,3};
          constant SI.Pressure henry_H[:] = {1.63e8,0} "Henry coefficient";
     constant Real henry_C[:]={2400,0}
  "constant to calculate temperature dependency";
     constant SI.Temperature henry_T_norm=298 "norm temperature";

  import
  ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices.ReferenceEnthalpy;


extends ThermalSeparation.Media.BaseMediumVapour(MMX={ 0.0440095,  0.01801528}, delta_hv_medium=true,nSubstance=2,V = {26.7,13.1},  eq_Tsonopoulos = {1,6}, sigma = {3.941,2.641}, epsilon_k = {195.2,809.1});
//   extends 
//     ThermalSeparation.Media.IdealGasMixtures.BaseClasses.BaseIdealGasMixture(delta_hv_medium = true,
//                               nS=2,nSubstance=2,V = {26.7,13.1},  eq_Tsonopoulos = {1,6}, sigma = {3.941,2.641}, epsilon_k = {195.2,809.1},
//     SpecificEnthalpy(start=if referenceChoice == ReferenceEnthalpy.ZeroAt0K then 
//                 3e5 else if referenceChoice == ReferenceEnthalpy.UserDefined then 
//                 h_offset else 0, nominal=1.0e5),
//     Density(start=2, nominal=10),
//     AbsolutePressure(start= 10e5, nominal=10e5),
//     Temperature(start=500, nominal=500),
//     fluidConstants={
//         ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.CO2,ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.H2O},
//                                                                           data=
//       {ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.CO2,ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.H2O}, R_const=data.R);


  redeclare replaceable model extends BaseProperties
    import
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasNasa;
          CalcSpecificEnthalpy calcSpecificEnthalpy(p=p,x=x,T=T);

  Real p_start(start=1.6e5);

    /*** Berechnung von eta ***/

  /***Berechnung der Sättigungstemperatur von Wasser beim Partialdruck des Wasserdampfes***/
protected
    Real pi "dimensionless pressure";
    Real[20] o "vector of auxiliary variables";
    SI.Temperature T_sat_water;

protected
   parameter Real conversionDebye = 3.33564e-30 "to convert from debye to C*m";
   SI.SpecificEnthalpy h_Jkg;
   SI.SpecificEnthalpy u_Jkg;
   // Real R;
  equation
  for i in 1:nSubstance loop
   X[i] = max(1e-5,min(1,x[i]*MMX[i]/MM));
  end for;

  p_start=state.p;
    //     assert(T >= 200 and T <= 6000, "
    // Temperature T (="   + String(T) + " K = 200 K) is not in the allowed range
    // 200 K <= T <= 6000 K
    // required from medium model \""   + mediumName + "\".");

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
    h_Jkg = X[2]*(4120*(T_sat_water-0) + 2382e3 + 1875*(T-T_sat_water))+ (X[1]*864)*(T-0);//h_TX(T, X);
     h = calcSpecificEnthalpy.h;
    // R = Modelica.Constants.R*X;
    u_Jkg = h_Jkg - Modelica.Constants.R*T;
    u = u_Jkg * MM;
    d = p/(Modelica.Constants.R*T);
    cp=2100*x[2]+919*x[1];
    eta=x[1]*0.00001847+x[2]*0.00001527;
    lambda=x[1]*0.02203+x[2]*0.02508;

    annotation (structurallyIncomplete);
  end BaseProperties;


  redeclare replaceable model extends CalcSpecificEnthalpy

  Real hX[nSubstance];

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
      hX[1] = 864*MMX[1] * (T-0);
     // hX[4] = 925*MMX[4] * (T- 0);
     // hX[1] = 1041*MMX[1]* (T-0);

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
  // for i in 3:nX loop
  //   h[i] = 0;
  // end for;

  end EvaporationEnthalpy;


  redeclare replaceable model extends FugacityCoefficient

  equation
  phi={1,1};

  annotation (Icon(graphics={
                        Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end FugacityCoefficient;


annotation (Icon(graphics={              Rectangle(
          extent={{-80,100},{100,-80}},
          lineColor={0,0,0},
          fillColor={215,230,240},
          fillPattern=FillPattern.Solid), Rectangle(
          extent={{-100,80},{80,-100}},
          lineColor={0,0,0},
          fillColor={240,240,240},
          fillPattern=FillPattern.Solid)}));
end CO2_H2O;
