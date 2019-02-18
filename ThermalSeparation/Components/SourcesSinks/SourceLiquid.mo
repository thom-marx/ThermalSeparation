within ThermalSeparation.Components.SourcesSinks;
model SourceLiquid
  extends Icons.Icons.LiquidSource;
replaceable package MediumLiquid =
    ThermalSeparation.Media.BaseMediumLiquid "medium to be used"                                    annotation(choicesAllMatching);
    MediumLiquid.BaseProperties medium(T0=T_ref, p=p, T=T_In_internal, x=x_in_internal,h=h);

  outer ThermalSeparation.SystemTS systemTS;
parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

parameter Boolean use_x_in = false
    "true if mole fraction variable is set by an input connector";

parameter ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption
    flowOption =                Enumerations.FlowOption.FlowVdot
    "determines which flow variable shall be specified"
    annotation(Evaluate=true,Dialog(group="Flow option"));

 final parameter Integer nS(min=2)=MediumLiquid.nSubstance;
  parameter SI.Temperature T=300 "temperature";

 // parameter SI.Concentration c[nS] = ones(nS);
  parameter Real Flow = 1
    "value of flow variable, for the applicable unit see flowOption"                       annotation(Dialog(group="Flow option"));
parameter SI.MoleFraction x[nS]= {0,1,0} "mole fraction";

parameter Boolean use_Flow = false
    "true if flow variable is set by an input connector"                                annotation(Dialog(group="Flow option"));
parameter Boolean useT_In = false
    "true if temperature variable is set by an input connector";
parameter Boolean fixed_pressure=false
    "true if pressure variable is set by an input connector";
parameter SI.Pressure p_fixed=1e5 "pressure" annotation (Dialog(enable=fixed_pressure));

Real c[nS];
SI.Pressure p;
ThermalSeparation.Units.MolarEnthalpy h(
                           start=-60e3);
SI.Density rho;
  ThermalSeparation.Interfaces.LiquidPortOut
                                    liquidPortOut(
                                         redeclare package Medium =
        MediumLiquid)
    annotation (Placement(transformation(extent={{114,0},{134,20}},   rotation=
            0), iconTransformation(extent={{94,-20},{134,20}})));
  Modelica.Blocks.Interfaces.RealInput Flow_In if use_Flow
    annotation (Placement(transformation(extent={{-120,20},{-80,60}}),
        iconTransformation(extent={{-120,20},{-80,60}})));

/*** monitoring ***/
Real X[nS] = x_in_internal[:].*MediumLiquid.MMX[:]/medium.MM;
//SI.MassFlowRate mdot = Vdot_In_internal * rho;
protected
 Modelica.Blocks.Interfaces.RealInput Flow_in_internal;
  Modelica.Blocks.Interfaces.RealInput T_In_internal;

public
  Modelica.Blocks.Interfaces.RealInput T_In if    useT_In
    annotation (Placement(transformation(extent={{-120,-20},{-80,20}}),
        iconTransformation(extent={{-120,-20},{-80,20}})));

protected
 Modelica.Blocks.Interfaces.RealInput x_in_internal[nS];
public
  Modelica.Blocks.Interfaces.RealInput x_in[nS] if use_x_in
    annotation (Placement(transformation(extent={{-120,-60},{-80,-20}},
                                                                      rotation=
            0), iconTransformation(extent={{-120,-60},{-80,-20}})));
equation
  if fixed_pressure then
    p=p_fixed;
  end if;
rho =  medium.d;
//medium.h = h;
c=x_in_internal/medium.MM*rho;

  if not useT_In then
  T_In_internal=T;
end if;
if not use_x_in then
  x_in_internal = x;
end if;
if not use_Flow then
     Flow_in_internal = Flow;
end if;

 if flowOption ==Enumerations.FlowOption.FlowVdot then
liquidPortOut.Ndot = -Flow_in_internal/medium.MM*rho;
 elseif flowOption ==Enumerations.FlowOption.FlowNdot then
   liquidPortOut.Ndot = -Flow_in_internal;
 else
    liquidPortOut.Ndot = -Flow_in_internal/medium.MM;
   end if;

 // liquidPortOut.T=T_In_internal;
 liquidPortOut.h_outflow = medium.h;
  liquidPortOut.x_outflow=x;

  //liquidPortOut.x = x_in_internal;
  liquidPortOut.p=p;

  connect(Flow_in_internal,Flow_In);
  connect(T_In_internal,T_In);
     for i in 1:nS loop
   connect(x_in_internal[i],x_in[i]);
     end for;

//        annotation (Diagram(graphics), Icon(graphics={Polygon(
//           points={{-40,72},{-80,-8},{0,-72},{80,-8},{40,72},{-40,72}},
//           lineColor={0,82,160},
//           smooth=Smooth.None,
//           fillColor={0,128,255},
//           fillPattern=FillPattern.Solid)}));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                   graphics), Documentation(info="<html>
<p>Liquid source. The user can choose whether to prescribe the mass flow rate, the molar flow rate or the volume flow rate. Therefore the enumeration parameter<i> flowOption</i> has to be set to the desired option. The value entered for the parameter <i>Flow</i> therefore corresponds either to mass flow rate in kg/s, molar flow rate in mol/s or volume flow rate in m3/s.</p>
<p>By default the flow rate, the temperature and the molar composition are prescribed but not the pressure. If the pressure is also to be prescribed, the parameter <i>p_fixed</i> has to be set to true.</p>
</html>"),
    Diagram(graphics));
end SourceLiquid;
