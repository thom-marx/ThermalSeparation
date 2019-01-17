within ThermalSeparation.Media;
package C2H5OH_Water_Vap "2 components: C2H5OH, H2O"
 constant SI.Temperature Tcrit_user[nS]= {513.9, 647.3};

       constant Integer henry[:] = {0};
     constant SI.Pressure henry_H[:] = {0,0} "Henry coefficient";
     constant Real henry_C[:]={0,0}
    "constant to calculate temperature dependency";
     constant SI.Temperature henry_T_norm=298 "norm temperature";

  import
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices.ReferenceEnthalpy;

        extends
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.BaseIdealGasMixture(delta_hv_medium = true,
                              nSubstance=2, V = {51.8, 13.1}, eq_Tsonopoulos = {4, 6},  sigma = {4.53, 2.641},epsilon_k = {362.6, 809.1},nS=2,MMX={0.046,0.018},
     SpecificEnthalpy(start=if referenceChoice == ReferenceEnthalpy.ZeroAt0K then
                 3e5 else if referenceChoice == ReferenceEnthalpy.UserDefined then
                 h_offset else 0, nominal=1.0e5),
     Density(start=2, nominal=10),
     AbsolutePressure(start= 10e5, nominal=10e5),
     Temperature(start=500, nominal=500),
     fluidConstants={ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.C2H5OH,
         ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.FluidData.H2O},
                                data=
       {ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.C2H5OH,
         ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasesData.H2O},
                                   R_const=data.R);
   redeclare replaceable model extends BaseProperties

    CalcSpecificEnthalpy calcSpecificEnthalpy(T0=T0,T=T, x=x, p=p);

       Real MM_min;

   EvaporationEnthalpy evapEnthalpy(p=p,T=T);

   equation
     for i in 1:nX loop
    X[i] = max(1e-5,min(1,x[i]*MMX[i]/MM));
   end for;
   d = p/(Modelica.Constants.R*T)*MM;
   MM_min=min(MMX);

   MM = max(MM_min,sum(x[j]*MMX[j] for j in 1:nS));
   //d=sum(c.*MMX);
   h = calcSpecificEnthalpy.h + sum(x.*evapEnthalpy.h);
     u = h - Modelica.Constants.R*T;

       R = data.R*X;

     annotation (Icon(graphics={
                           Rectangle(
               extent={{-100,100},{100,-100}},
               lineColor={0,0,255},
               fillColor={255,255,255},
               fillPattern=FillPattern.Solid)}));
   end BaseProperties;

  model Testmodell
    parameter Integer n=1;
    parameter Integer nS=2;
  parameter SI.Temperature T = 373.15;
  parameter SI.Pressure p = 100000;
  parameter SI.MoleFraction x[nS] = {0.3,0.7};//,0.05}};
  parameter SI.Concentration c[nS]={1,1};

  BaseProperties BP(T=T,p=p,x=x,c=c, x_star=x);
  end Testmodell;

  redeclare replaceable model extends EvaporationEnthalpy

  protected
        constant ThermalSeparation.Units.MolarEnthalpy h_v_H2O =  2382.2e3*MMX[2];// spezifische Verdampfungsenthalpie bei 50°C
           ThermalSeparation.Units.MolarEnthalpy h_v_EtOH;
  equation
      /*** enthalpy ***/
    h[1] = h_v_EtOH;
    h[2] = h_v_H2O;

    h_v_EtOH = (50.43*exp(0.4475*(T/Tcrit_user[1]))*(1-(T/Tcrit_user[1]))^0.4989)*1000;  // J/mol

  end EvaporationEnthalpy;

  redeclare replaceable model extends CalcSpecificEnthalpy

    /***Berechnung der Sättigungstemperatur von Wasser beim Partialdruck des Wasserdampfes***/

   constant SI.MolarHeatCapacity cp_v_EtOH = 1765*MMX[1]; // J/kgK
    constant SI.MolarHeatCapacity cp_v_H2O = 2027*MMX[2];
    constant SI.MolarHeatCapacity cp_l_EtOH = 2443.48*MMX[1];
    constant SI.MolarHeatCapacity cp_l_H2O = 4187*MMX[2];
     // constant ThermalSeparation.Units.MolarEnthalpy h_v_H2O =  2382.2e3*MMX[2];// spezifische Verdampfungsenthalpie bei 50°C

      SI.MolarMass MM;

   //  ThermalSeparation.Units.MolarEnthalpy h_v_EtOH;

  protected
         ThermalSeparation.Units.MolarEnthalpy h_l_Jkg;
     ThermalSeparation.Units.MolarEnthalpy h_v_Jkg;
     ThermalSeparation.Units.MolarEnthalpy h_l_EtOH_Jkg;
     ThermalSeparation.Units.MolarEnthalpy h_l_H2O_Jkg;
     ThermalSeparation.Units.MolarEnthalpy h_v_EtOH_Jkg;
     ThermalSeparation.Units.MolarEnthalpy h_v_H2O_Jkg;

  equation
    MM = sum(x[j]*MMX[j] for j in 1:nS);

   // h_v_EtOH = (50.43*exp(0.4475*(T/Tcrit_user[1]))*(1-(T/Tcrit_user[1]))^0.4989)*1000;  // J/mol

    h_l_EtOH_Jkg=cp_l_EtOH*(351.5-T0); //T_boil ethanol constant 351.5 K
  h_l_H2O_Jkg=cp_l_H2O*(373.15-T0);  // T boil water constant 373.15 K
  h_v_EtOH_Jkg=cp_v_EtOH*(T-351.5);
  h_v_H2O_Jkg=cp_v_H2O*(T-373.15);
  h_l_Jkg=x[1]*h_l_EtOH_Jkg+x[2]*h_l_H2O_Jkg;
  h_v_Jkg=x[1]*h_v_EtOH_Jkg+x[2]*h_v_H2O_Jkg;

  h =h_l_Jkg +  h_v_Jkg; ///

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;

  annotation (Documentation(info="<html>
  
</html>"));
end C2H5OH_Water_Vap;
