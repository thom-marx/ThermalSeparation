within ThermalSeparation.Components.SourcesSinks;
model CombGas_x
     extends Icons.Icons.LiquidMixer;
  final parameter Integer nS=Medium.nSubstance;

  outer ThermalSeparation.SystemTS systemTS;
parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

  replaceable package Medium =
      ThermalSeparation.Media.BaseMediumVapour constrainedby
    ThermalSeparation.Media.BaseMediumVapour "medium to be used"                                                         annotation(choicesAllMatching);
    Medium.BaseProperties medium(T0=T_ref,p=1.5e5, T=T, x=x, x_star=x, c=c);
SI.Temperature T;
parameter Real mass = 10 "mass";
parameter SI.Volume V = 0.01 "volume";
parameter Real[nS] x_start= {0,1,0} "initial value for mole fraction";
Real MM_start(start=0.018) = medium.MM "start value for molar mass";
parameter SI.Temperature T_start=50+273 "initial value for Temperature";
//Real c[Medium.nSubstance]= x/medium.MM*medium.d;
  Interfaces.GasPortIn    gasPortIn2(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{-100,-50},{-80,-30}},rotation=
           0), iconTransformation(extent={{-120,-70},{-80,-30}})));
  Interfaces.GasPortIn    gasPortIn1(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{-100,50},{-80,70}},rotation=0),
        iconTransformation(extent={{-120,30},{-80,70}})));
  Interfaces.GasPortOut    gasPortOut(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{100,0},{120,20}}, rotation=0),
        iconTransformation(extent={{80,-20},{120,20}})));

Real c[Medium.nSubstance];
Real x[Medium.nSubstance](start=x_start) = c*medium.MM/medium.d;
equation
  gasPortIn1.h_outflow = medium.h;
  gasPortIn2.h_outflow = medium.h;
  gasPortOut.h_outflow = medium.h;
  gasPortIn1.x_outflow = x;
  gasPortIn2.x_outflow = x;
  gasPortOut.x_outflow = x;

  /*** energy balance ***/
gasPortOut.Ndot * gasPortOut.h_outflow  +  gasPortIn1.Ndot* inStream(gasPortIn1.h_outflow) + gasPortIn2.Ndot* inStream(gasPortIn2.h_outflow) = V* der(medium.u* sum(c));
/*** mole balance ***/
gasPortOut.Ndot * gasPortOut.x_outflow +  inStream(gasPortIn1.x_outflow) * gasPortIn1.Ndot + inStream(gasPortIn2.x_outflow) * gasPortIn2.Ndot =V*der(c);

//gasPortOut.Vdot * gasPortOut.x * sum(gasPortOut.c) +  gasPortIn1.x * gasPortIn1.Vdot* sum(gasPortIn1.c) + gasPortIn2.x * gasPortIn2.Vdot* sum(gasPortIn2.c) = zeros(nS);

gasPortOut.p = gasPortIn1.p;
gasPortOut.p = gasPortIn2.p;

gasPortIn1.Ndot + gasPortIn2.Ndot + gasPortOut.Ndot = V*der(medium.d/medium.MM);

initial equation
  T = T_start;
 c = x_start/medium.MM*medium.d;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end CombGas_x;
