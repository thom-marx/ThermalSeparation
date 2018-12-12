within ThermalSeparation.Components.SourcesSinks;
model SourceGas "source for gas"
 extends Icons.Icons.GasSource;
replaceable package Medium =
   ThermalSeparation.Media.BaseMediumVapour "medium to be used"                                                     annotation(choicesAllMatching);
Medium.BaseProperties medium(p=gasPortOut.p,x=x_in_internal, T=T_in_internal, c=c, x_star=x_in_internal);
parameter Boolean use_T_in= false
    "true if temperature variable is set by an input connector";
parameter Boolean use_x_in = false
    "true if mole fraction variable is set by an input connector";
parameter Boolean use_Flow = true
    "true if flow variable is set by an input connector";
 // parameter Integer nS(min=1)=2;
  parameter SI.Temperature T = 300 "temperature";

parameter ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption
    flowOption =                Enumerations.FlowOption.FlowVdot
    "determines which flow variable shall be specified"
    annotation(Evaluate=true);

  parameter Real Flow= 0.5
    "value of flow variable, for the applicable unit see flowOption"                        annotation (Dialog(enable= (flowOption == Enumerations.FlowOption.FlowVdot)));
 // parameter SI.MolarFlowRate Ndot = 1 annotation (Dialog(enable= ( flowOption == FlowOption.FlowNdot)));
  parameter SI.MoleFraction x[Medium.nSubstance]= {1,0} "mole fraction";
  ThermalSeparation.Interfaces.GasPortOut
                                 gasPortOut(
                                   redeclare package Medium=Medium)       annotation (Placement(transformation(
          extent={{114,0},{134,20}}, rotation=0), iconTransformation(extent={{
            94,-20},{134,20}})));

Real rho;
Real c[Medium.nSubstance];
  ThermalSeparation.Units.MolarEnthalpy h(
                             start=430e3);
  Modelica.Blocks.Interfaces.RealInput Flow_in if use_Flow
    annotation (Placement(transformation(extent={{-120,20},{-80,60}}, rotation=
            0), iconTransformation(extent={{-120,20},{-80,60}})));
  Modelica.Blocks.Interfaces.RealInput T_in if use_T_in
    annotation (Placement(transformation(extent={{-120,-20},{-80,20}},rotation=
            0), iconTransformation(extent={{-120,-20},{-80,20}})));
  Modelica.Blocks.Interfaces.RealInput x_in[Medium.nSubstance] if use_x_in
    annotation (Placement(transformation(extent={{-120,-60},{-80,-20}},
                                                                      rotation=
            0), iconTransformation(extent={{-120,-60},{-80,-20}})));

SI.MolarMass MM(start=0.02) = medium.MM;
protected
  Modelica.Blocks.Interfaces.RealInput T_in_internal;
  Modelica.Blocks.Interfaces.RealInput x_in_internal[Medium.nSubstance];
  Modelica.Blocks.Interfaces.RealInput Flow_in_internal;

equation
  if not use_T_in then
    T_in_internal = T;
  end if;
  if not use_x_in then
    x_in_internal = x;
  end if;
  if not use_Flow then
    Flow_in_internal = Flow;
  end if;

medium.h=h;

rho = medium.d;
c=x_in_internal/medium.MM*rho;

 //gasPortOut.T=T_in_internal;
 gasPortOut.h_outflow = medium.h;
 gasPortOut.x_outflow=x;
 if flowOption ==Enumerations.FlowOption.FlowVdot then
gasPortOut.Ndot = -Flow_in_internal/MM*rho;
 elseif flowOption ==Enumerations.FlowOption.FlowNdot then
   gasPortOut.Ndot = -Flow_in_internal;
 else
    gasPortOut.Ndot = -Flow_in_internal/MM;
   end if;
//gasPortOut.x = x_in_internal;
//gasPortOut.p_medium = gasPortOut.p;
   connect(T_in_internal,T_in);
   for i in 1:Medium.nSubstance loop
   connect(x_in_internal[i],x_in[i]);
   end for;
   connect(Flow_in_internal,Flow_in);

  annotation (Diagram(graphics), Icon(coordinateSystem(preserveAspectRatio=
            false, extent={{-100,-100},{100,100}}), graphics));
end SourceGas;
