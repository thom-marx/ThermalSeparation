within ThermalSeparation.Components.SourcesSinks;
model SplitLiquid_1p "splitter, defining outlet pressure at only one port"
   extends Icons.Icons.LiquidSplitter;
   replaceable package Medium = Media.BaseMediumLiquid annotation(choicesAllMatching);
  final parameter Integer nS(min=1)=Medium.nSubstance;

  ThermalSeparation.Interfaces.LiquidPortIn
                          liquidPortIn(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{100,0},{120,20}},   rotation=
            0), iconTransformation(extent={{80,-20},{120,20}})));
  ThermalSeparation.Interfaces.LiquidPortOut
                           liquidPortOut1(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{-100,-50},{-80,-30}},
                                                                   rotation=0),
        iconTransformation(extent={{-120,-70},{-80,-30}})));
  ThermalSeparation.Interfaces.LiquidPortOut
                           liquidPortOut2(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{-100,50},{-80,70}},rotation=0),
        iconTransformation(extent={{-120,30},{-80,70}})));

  parameter SourcesSinks.Enumerations.SplitOption splitOption=Enumerations.SplitOption.out1_over_in
    "split options"                                                                      annotation(Dialog(Evaluate=true));

  Modelica.Blocks.Interfaces.RealInput split
    annotation (Placement(transformation(extent={{-94,52},{-54,92}}),
        iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={0,60})));

//Real x[Medium.nSubstance](start=x_start);//= c*medium.MM/medium.d;
//Real h;
parameter SI.Volume V = 0.01;

equation
//liquidPortIn.p = 1e5;

  liquidPortOut1.h_outflow = inStream(liquidPortIn.h_outflow);
  liquidPortOut2.h_outflow = inStream(liquidPortIn.h_outflow);
  liquidPortIn.h_outflow = 4;
  liquidPortOut1.x_outflow = inStream(liquidPortIn.x_outflow);
  liquidPortOut2.x_outflow = inStream(liquidPortIn.x_outflow);
  liquidPortIn.x_outflow = fill(1/Medium.nSubstance,Medium.nSubstance);

liquidPortOut2.p = liquidPortIn.p;

liquidPortOut1.Ndot + liquidPortOut2.Ndot + liquidPortIn.Ndot = 0;

 if splitOption==Enumerations.SplitOption.out1_over_in then
liquidPortOut1.Ndot + split *liquidPortIn.Ndot = 0;
 elseif splitOption==Enumerations.SplitOption.out2_over_in then
  liquidPortOut2.Ndot + split *liquidPortIn.Ndot = 0;
 elseif splitOption==Enumerations.SplitOption.out2_over_out1 then
   split*liquidPortOut1.Ndot - liquidPortOut2.Ndot=0;
 elseif splitOption==Enumerations.SplitOption.Nout1_over_Nin then
liquidPortOut1.Ndot + split *liquidPortIn.Ndot  = 0;
 end if;
// initial equation
//     T = 50+273;
//  x = x_start;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end SplitLiquid_1p;
