within ThermalSeparation.Components;
package Condenser 
  extends Icons.Library.Red;

  model FlashCondenser_CO2_H2O "for CO2/H2O: cools down and completely separates streams without considering energy removal"

  // Check luft nicht, weil die Vektorgre der Outputs hier fest vorgegeben ist und nicht mit den Defaults der Connectoren bereinstimmt. Das macht aber auch erstmal nichts, weil die Medienmodell hier nicht replaceable sind.

      extends ThermalSeparation.Icons.Color.Condenser;

  // parameter Integer nSV = 2;
  // parameter Integer nSL = 3;
  parameter Integer onlyVap = 2 "component only in vapour stream";
  parameter Integer onlyLiq = 1 "component only in liquid stream";

  parameter Modelica.SIunits.Temperature T_out=273.15 + 40 "temperature at outlet";

  Modelica.SIunits.Density rho_v_in=mediumVapourIn.d;
  Modelica.SIunits.Density rho_v=mediumVapour.d;
  Modelica.SIunits.Density rho_l=mediumLiquid.d;

  Modelica.SIunits.MolarMass MM_v_in(start=0.028)=mediumVapourIn.MM;
  Modelica.SIunits.MolarMass MM_v(start=0.044)=mediumVapour.MM;
  Modelica.SIunits.MolarMass MM_l=mediumLiquid.MM;

     package MediumVapour =
        ThermalSeparation.Media.H2O_CO2_Vap;

      MediumVapour.BaseProperties mediumVapourIn(
      T0=293.15,
      p=gasPortIn.p,
      T=T_v_in,
      x=inStream(gasPortIn.x_outflow),
      c=inStream(gasPortIn.x_outflow)*rho_v_in/MM_v_in, x_star=inStream(gasPortIn.x_outflow));

      MediumVapour.BaseProperties mediumVapour(
      T0=293.15,
      p=gasPortOut.p,
      T=T_out,
      x=gasPortOut.x_outflow,
      c=gasPortOut.x_outflow*rho_v/MM_v, x_star=gasPortOut.x_outflow);

     package MediumLiquid =
      ThermalSeparation.Media.H2O_CO2_MEA_Liq;

     MediumLiquid.BaseProperties mediumLiquid(
      T0=293.15,
      p=gasPortOut.p,
      T=T_out,
      x=liquidPortOut.x_outflow,
      h=gasPortIn.h_outflow);

      Modelica.SIunits.HeatFlowRate Q;

     //  Real checkbal;
  Real h_v_in = inStream(gasPortIn.h_outflow);
  Real T_v_in;
    ThermalSeparation.Interfaces.GasPortIn gasPortIn(redeclare package Medium =
          MediumVapour) annotation (Placement(transformation(extent={{-100,40},{-80,
              60}}), iconTransformation(extent={{-112,-20},{-72,20}})));
    ThermalSeparation.Interfaces.GasPortOut gasPortOut(redeclare package Medium =
          MediumVapour) annotation (Placement(transformation(extent={{-100,60},{-80,
              80}}), iconTransformation(extent={{-24,66},{16,106}})));
    ThermalSeparation.Interfaces.LiquidPortOut liquidPortOut(redeclare package Medium =
                 MediumLiquid) annotation (Placement(transformation(extent={{40,-100},
              {60,-80}}), iconTransformation(extent={{-26,-106},{16,-66}})));

    ThermalSeparation.Interfaces.HeatPort heatPort(Qdot=-Q) annotation (Placement(
          transformation(extent={{-100,-100},{-80,-80}}, rotation=0),
          iconTransformation(extent={{60,-20},{100,20}})));
  equation
       h_v_in = mediumVapourIn.h;
       //gasPortIn.h_outflow=mediumLiquid.h;
       gasPortIn.x_outflow={1,0};

     gasPortOut.p = gasPortIn.p;
    gasPortOut.x_outflow = {0,1};
    -gasPortOut.Ndot * gasPortOut.x_outflow[onlyVap] = inStream(gasPortIn.x_outflow[onlyVap]) * gasPortIn.Ndot;
    gasPortOut.h_outflow = mediumVapour.h;

    liquidPortOut.x_outflow = {1,0,0};
    if time>=5 then
    -liquidPortOut.Ndot * liquidPortOut.x_outflow[onlyLiq] = inStream(gasPortIn.x_outflow[onlyLiq]) * gasPortIn.Ndot;
    else
     liquidPortOut.Ndot=0;
    end if;
    liquidPortOut.h_outflow = mediumLiquid.h;

   Q = h_v_in * gasPortIn.Ndot + mediumLiquid.h * liquidPortOut.Ndot + mediumVapour.h * gasPortOut.Ndot;
  //  checkbal = gasPortIn.Vdot * sum(gasPortIn.c[:]) + gasPortOut.Vdot * sum(gasPortOut.c[:]) + liquidPortOut.Vdot * sum(liquidPortOut.c[:]);

   annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                    graphics), Diagram(graphics));
  end FlashCondenser_CO2_H2O;
end Condenser;
