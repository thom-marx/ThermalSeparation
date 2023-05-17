within ThermalSeparation.Components.SourcesSinks;
model CombLiquid_x2
     extends Icons.Color.LiquidMixer;
  final parameter Integer nS=Medium.nSubstance;

  outer ThermalSeparation.SystemTS systemTS;
parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

  replaceable package Medium =
      ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O             constrainedby ThermalSeparation.Media.BaseMediumLiquid
                                             "medium to be used"                                                         annotation(choicesAllMatching);
    Medium.BaseProperties medium(T0=T_ref,p=liquidPortOut.p, T=T, x=x,h=liquidPortOut.h_outflow);
    SI.Temperature T;

parameter Real mass = 10 "mass";
parameter SI.Volume V = 0.01 "volume";
parameter Real[nS] x_start= {0,1,0} "initial value for mole fraction";
Real MM_start(start=0.018) = medium.MM "start value for molar mass";
parameter SI.Temperature T_start=50+273 "initial value for Temperature";

  ThermalSeparation.Interfaces.LiquidPortIn
                          liquidPortIn2(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{-100,-50},{-80,-30}},rotation=
           0), iconTransformation(extent={{-120,-70},{-80,-30}})));
  ThermalSeparation.Interfaces.LiquidPortIn
                          liquidPortIn1(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{-100,50},{-80,70}},rotation=0),
        iconTransformation(extent={{-120,30},{-80,70}})));
  ThermalSeparation.Interfaces.LiquidPortOut
                           liquidPortOut(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{100,0},{120,20}}, rotation=0),
        iconTransformation(extent={{80,-20},{120,20}})));

Real c[Medium.nSubstance];
Real x[Medium.nSubstance](start=x_start) = c*medium.MM/medium.d;
equation
  liquidPortIn1.h_outflow = medium.h;
  liquidPortIn2.h_outflow = medium.h;
  //liquidPortOut.h_outflow = medium.h;
  liquidPortIn1.x_outflow = x;
  liquidPortIn2.x_outflow = x;
  liquidPortOut.x_outflow = x;

  /*** energy balance ***/
liquidPortOut.Ndot * liquidPortOut.h_outflow  +  liquidPortIn1.Ndot* inStream(liquidPortIn1.h_outflow) + liquidPortIn2.Ndot* inStream(liquidPortIn2.h_outflow) = V* der(medium.u* sum(c));
/*** mole balance ***/
liquidPortOut.Ndot * liquidPortOut.x_outflow +  inStream(liquidPortIn1.x_outflow) * liquidPortIn1.Ndot + inStream(liquidPortIn2.x_outflow) * liquidPortIn2.Ndot =V*der(c);

//liquidPortOut.Vdot * liquidPortOut.x * sum(liquidPortOut.c) +  liquidPortIn1.x * liquidPortIn1.Vdot* sum(liquidPortIn1.c) + liquidPortIn2.x * liquidPortIn2.Vdot* sum(liquidPortIn2.c) = zeros(nS);

liquidPortOut.p = liquidPortIn1.p;
liquidPortOut.p = liquidPortIn2.p;

liquidPortIn1.Ndot + liquidPortIn2.Ndot + liquidPortOut.Ndot = V*der(medium.d/medium.MM);

initial equation
  T = T_start;
 c = x_start/medium.MM*medium.d;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end CombLiquid_x2;
