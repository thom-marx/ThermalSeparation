within ThermalSeparation.Media.IdealGasMixtures;
package Propane_Pentane "Propane, Pentane"

/** henry's law is not applicable for this medium **/

     constant Integer henry[:] = {0,0};
          constant SI.Pressure henry_H[:] = {0,0} "Henry coefficient";
     constant Real henry_C[:]={0,0}
  "constant to calculate temperature dependency";
     constant SI.Temperature henry_T_norm=298 "norm temperature";

  import
  ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices.ReferenceEnthalpy;


extends ThermalSeparation.Media.BaseMediumVapour(MMX={ 0.0441,  0.07215}, delta_hv_medium=true,nSubstance=2,V = {66.18,107.22},  eq_Tsonopoulos = {1,1}, sigma = {0.016,0.016}, epsilon_k = {206,269});

// ref. for epsilon_k (lennard jones force constant): http://onlinelibrary.wiley.com/doi/10.1002/aic.690080320/pdf


  redeclare replaceable model extends BaseProperties
    import
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasNasa;
          CalcSpecificEnthalpy calcSpecificEnthalpy(p=p,x=x,T=T);

  Real p_start(start=1.6e5);

protected
    SI.Temperature T_sat[nSubstance]={231.1,309.2};

protected
   parameter Real conversionDebye = 3.33564e-30 "to convert from debye to C*m";
    Real R;
  equation
  for i in 1:nSubstance loop
   X[i] = max(1e-5,min(1,x[i]*MMX[i]/MM));
  end for;

  p_start=state.p;

    //bei der Gleichung MM=molarMass(state): Index-Problem!
   MM = max(1e-6,sum(x.*MMX));//molarMass(state);
    //h_Jkg = X[2]*(4120*(T_sat-0) + 2382e3 + 1875*(T-T_sat))+ (X[1]*864)*(T-0);//h_TX(T, X);
     h = calcSpecificEnthalpy.h;
     R = Modelica.Constants.R/MM;
     //R = data.R_s*X;
    //u_Jkg = h_Jkg - Modelica.Constants.R*T;
    u = h-R*T*MM;
    //u = u_Jkg * MM;
    d = p/(R*T);
    cp=1549*x[1]+2456*x[2];
    eta=7.7e-6;
    lambda=x[1]*0.018+x[2]*0.0144;

    annotation (structurallyIncomplete);
  end BaseProperties;

  redeclare replaceable model extends CalcSpecificEnthalpy

  Real hX[nSubstance];

protected
   SI.Temperature T_sat[nSubstance]={231.1,309.2};

  equation
      h = sum(x.*hX);
  //     hX[1] = (2456*(T_sat[1]-0) + 426e3 + 1549*(T-T_sat[1]))*MMX[1];
  //     hX[2] = (2206*(T_sat[2]-0) + 358e3 + 1599*(T-T_sat[2]))*MMX[2];
      hX[1] = (2456*(T_sat[1]-1406) + 426e3 + 1549*(T-T_sat[1]))*MMX[1];
      hX[2] = (2206*(T_sat[2]-1406) + 358e3 + 1599*(T-T_sat[2]))*MMX[2];

  //reference temperature T0=1406 has been modified so that h matches the value of FP

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;

  redeclare replaceable model extends EvaporationEnthalpy

   // SI.SpecificEnthalpy r_water=-2462.5*T + 3177.8e3;
  equation
      /*** enthalpy ***/
      h[1]=426e3*MMX[1];
    h[2] = 358e3*MMX[2];

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
end Propane_Pentane;
