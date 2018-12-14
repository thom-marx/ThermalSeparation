within ThermalSeparation.Components.Pumps;
model idealPumpControlledVdot
 extends Icons.Icons.Pumps;
  outer ThermalSeparation.SystemTS systemTS;
parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

replaceable package MediumLiquid =
ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O     constrainedby
    Media.BaseMediumLiquid                                                                            annotation(choicesAllMatching);
  MediumLiquid.BaseProperties mediumLiquid(T0=T_ref, p=p, T=T_l, x=x_l,h=h_l);
    MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref, p=p_in, T=T_l_in, x=x_l_in,h=h_l_in);

     ThermalSeparation.Interfaces.LiquidPortIn
                             liquidIn(redeclare package Medium=MediumLiquid)
                                        annotation (Placement(transformation(
          extent={{2,-94},{22,-74}},   rotation=0), iconTransformation(extent={{-22,
            -114},{22,-74}})));

   ThermalSeparation.Interfaces.LiquidPortOut
                            liquidOut(redeclare package Medium=MediumLiquid)
                                          annotation (Placement(transformation(
          extent={{0,102},{20,122}},
                                   rotation=0), iconTransformation(extent={{-20,82},
            {20,122}})));

parameter SI.Time T = 10 "Delay time for outlet volume flow rate V_flow_Out";
parameter SI.TemperatureDifference dT = 0.5 "Temperature increase T_Out - T_In";
parameter SI.VolumeFlowRate V_flow_start=0.000305
    "Fixed initial value for outlet volume flow rate V_flow_Out";

  final parameter Integer nSL = MediumLiquid.nSubstance;

/*** Medium properties ***/
MediumLiquid.Density rho_l;
SI.Concentration c_l[nSL](each start = 500, each nominal = 1e4);
SI.MoleFraction x_l[nSL];
  ThermalSeparation.Units.MolarEnthalpy h_l;
  ThermalSeparation.Units.MolarEnthalpy u_l;
  SI.MolarMass MM_l;

SI.MoleFraction x_l_in[nSL];
  ThermalSeparation.Units.MolarEnthalpy h_l_in;//=mediumLiquidIn.h;
  MediumLiquid.Temperature T_l_in;
  MediumLiquid.Temperature T_l(start=350);
  SI.VolumeFlowRate Vdot_l(start=1e-3, nominal = 1e-3)= liquidOut.Ndot *MM_l/rho_l;
  SI.VolumeFlowRate Vdot_l_in(start=1e-3, nominal = 1e-3)=liquidIn.Ndot *mediumLiquidIn.MM/mediumLiquidIn.d;

  MediumLiquid.AbsolutePressure p_out;
  MediumLiquid.AbsolutePressure p_in;
  MediumLiquid.AbsolutePressure p(start=2e5);
  // SI.VolumeFlowRate V_flow_delayed(stateSelect=StateSelect.default);

/*** electric power consumption ***/
parameter Boolean powerCons = false
    "true if power consumption shall be calculated using nominal values";
parameter SI.VolumeFlowRate Vdot_nom(min=1e-10) = 1
    "nominal volume flow rate for power consumption" annotation(Dialog(enable=powerCons));
parameter SI.Power P_nom= 0 "nominal power consumption" annotation(Dialog(enable=powerCons));
SI.Power P = P_nom * Vdot_l_in^3/ Vdot_nom^3 "power consumption";

/*** monitoring***/
SI.MassFlowRate mdot=-liquidOut.Ndot*MM_l;

  Modelica.Blocks.Interfaces.RealInput VdotInput
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}}),
        iconTransformation(extent={{-140,-20},{-100,20}})));
initial equation
//V_flow_delayed = V_flow_start;

equation
liquidIn.x_outflow = x_l;
liquidIn.h_outflow = h_l;
  //der(V_flow_delayed) = (Vdot_l-V_flow_delayed)/T;
   //downstream
   liquidIn.p = p_in;
   liquidOut.p = p_out;

   inStream(liquidIn.h_outflow) = h_l_in;
   liquidOut.h_outflow = h_l;
   inStream(liquidIn.x_outflow) = x_l_in;
   liquidOut.x_outflow = x_l;

  //h_l =mediumLiquid.h;
  u_l =mediumLiquid.u;
   rho_l = mediumLiquid.d;
   MM_l = mediumLiquid.MM;

  /*** correlation between mole fraction x and concentration c ***/
  for i in 1:nSL loop
    x_l[i] = c_l[i] * MM_l /rho_l;
  end for;

    /***energy balance ***/
 //0  = Vdot_l_in * sum(c_l_in[:])*h_l_in - Vdot_l*sum(c_l[:])*h_l;
  T_l = T_l_in+dT;

    /*** mass balance ***/
   // for i in 1:nSL loop
          // A*H* der(c_l[i]) = Vdot_l_in*c_l_in[i]  - Vdot_l * c_l[i];
        liquidIn.Ndot + liquidOut.Ndot = 0;
          Vdot_l = -VdotInput;
   // end for;

   x_l = x_l_in;

  p_out=p;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                   graphics));
end idealPumpControlledVdot;
