within ThermalSeparation.Components.SourcesSinks;
model SplitGas_1p "splitter, defining outlet pressure at only one port"
   extends Icons.Icons.LiquidSplitter;
   replaceable package Medium = Media.BaseMediumVapour annotation(choicesAllMatching);
  final parameter Integer nS(min=1)=Medium.nSubstance;

  Interfaces.GasPortIn gasPortIn(redeclare package Medium = Medium) annotation (
     Placement(transformation(extent={{100,0},{120,20}}, rotation=0),
        iconTransformation(extent={{80,-20},{120,20}})));
  Interfaces.GasPortOut gasPortOut1(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-100,-50},{-80,-30}},
          rotation=0), iconTransformation(extent={{-120,-70},{-80,-30}})));
  Interfaces.GasPortOut gasPortOut2(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-100,50},{-80,70}}, rotation=0),
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

  gasPortOut1.h_outflow = inStream(gasPortIn.h_outflow);
  gasPortOut2.h_outflow = inStream(gasPortIn.h_outflow);
  gasPortIn.h_outflow = 4;
  gasPortOut1.x_outflow = inStream(gasPortIn.x_outflow);
  gasPortOut2.x_outflow = inStream(gasPortIn.x_outflow);
  gasPortIn.x_outflow = fill(1/Medium.nSubstance,Medium.nSubstance);

gasPortOut2.p = gasPortIn.p;

gasPortOut1.Ndot + gasPortOut2.Ndot + gasPortIn.Ndot = 0;

 if splitOption==Enumerations.SplitOption.out1_over_in then
gasPortOut1.Ndot + split *gasPortIn.Ndot = 0;
 elseif splitOption==Enumerations.SplitOption.out2_over_in then
  gasPortOut2.Ndot + split *gasPortIn.Ndot = 0;
 elseif splitOption==Enumerations.SplitOption.out2_over_out1 then
   split*gasPortOut1.Ndot - gasPortOut2.Ndot=0;
 elseif splitOption==Enumerations.SplitOption.Nout1_over_Nin then
gasPortOut1.Ndot + split *gasPortIn.Ndot  = 0;
 end if;
// initial equation
//     T = 50+273;
//  x = x_start;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end SplitGas_1p;
