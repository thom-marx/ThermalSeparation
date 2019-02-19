within ThermalSeparation.Components.SourcesSinks;
model SinkLiquid_medium
   extends Icons.Color.LiquidSink;
replaceable package MediumLiquid =
    ThermalSeparation.Media.BaseMediumLiquid "medium to be used"                                    annotation(choicesAllMatching);
//MediumLiquid.BaseProperties medium(T0=T_ref, p=p, T=T, x=inStream(liquidPortIn.x_outflow)); // parameter Integer nS(min=1);
  MediumLiquid.BaseProperties medium(T0=T_ref, p=p, T=T, x={0.5,0.5},h=h_liq);
  //parameter SI.Length h = 5;
  parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));
 outer ThermalSeparation.SystemTS systemTS;
 SI.Temperature T;
 SI.MolarEnthalpy h_liq;
  ThermalSeparation.Interfaces.LiquidPortIn
                                   liquidPortIn(
                                       redeclare package Medium=MediumLiquid) annotation (Placement(transformation(
          extent={{-96,0},{-76,20}}, rotation=0), iconTransformation(extent={{-116,
            -20},{-76,20}})));
parameter Boolean use_p = true "fixed pressure";
parameter SI.Pressure p = 1.3e5 "pressure" annotation(Dialog(enable=use_p));

parameter Boolean use_Vdot = false "fixed volume flow rate";

  parameter Boolean showVdot = false "real output showing Vdot";
   Modelica.Blocks.Interfaces.RealOutput Vdot = liquidPortIn.Ndot if showVdot annotation (Placement(
         transformation(extent={{66,-6},{100,28}}), iconTransformation(extent={{-17,-17},
            {17,17}},
        rotation=180,
        origin={83,-83})));

   Modelica.Blocks.Interfaces.RealInput Vdot_set if use_Vdot
     annotation (Placement(transformation(extent={{-100,20},{-60,60}}, rotation=
             0), iconTransformation(extent={{-15,-15},{15,15}},
        rotation=180,
        origin={85,85})));

protected
   Modelica.Blocks.Interfaces.RealInput Vdot_set_internal;

equation
  medium.h = inStream(liquidPortIn.h_outflow);
if use_p then
  liquidPortIn.p = p;
end if;

 if use_Vdot then
   liquidPortIn.Ndot = Vdot_set_internal;
 else
   Vdot_set_internal=1;
 end if;

  liquidPortIn.x_outflow = fill(1/MediumLiquid.nSubstance,MediumLiquid.nSubstance);
  liquidPortIn.h_outflow = 1;

      connect(Vdot_set_internal,Vdot_set);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                   graphics));
end SinkLiquid_medium;
