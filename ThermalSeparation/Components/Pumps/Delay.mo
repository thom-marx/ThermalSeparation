within ThermalSeparation.Components.Pumps;
model Delay
replaceable package Medium = Media.BaseMediumLiquid annotation(choicesAllMatching);
    outer ThermalSeparation.SystemTS systemTS;
   final parameter Integer nSL = Medium.nSubstance "number of substances";

   Medium.BaseProperties medium(p=liquidIn.p, x= inStream(liquidIn.x_outflow),h = liquidOut.h_outflow, T=T_med);

     ThermalSeparation.Interfaces.LiquidPortIn
                             liquidIn(redeclare package Medium=Medium)
                                        annotation (Placement(transformation(
          extent={{-100,-14},{-80,6}}, rotation=0), iconTransformation(extent={{-100,
            -14},{-80,6}})));

   ThermalSeparation.Interfaces.LiquidPortOut
                            liquidOut(redeclare package Medium=Medium)
                                          annotation (Placement(transformation(
          extent={{80,-8},{100,12}},
                                   rotation=0), iconTransformation(extent={{80,-8},
            {100,12}})));

parameter SI.Time T = 10 "Delay time for outlet volume flow rate V_flow_Out";
parameter SI.VolumeFlowRate V_flow_start=0.000305
    "Fixed initial value for outlet volume flow rate V_flow_Out";

  SI.MolarFlowRate Ndot_l(start=1e-3, nominal = 1e-3);

  SI.MolarFlowRate N_flow_delayed(stateSelect=StateSelect.default);
  SI.Temperature T_med;

initial equation
N_flow_delayed = V_flow_start * medium.d/medium.MM;

equation
  //medium.h = liquidOut.h_outflow;
   liquidIn.h_outflow = inStream(liquidOut.h_outflow);
   liquidIn.x_outflow = inStream(liquidOut.x_outflow);

  der(N_flow_delayed) = (Ndot_l-N_flow_delayed)/T;
   //downstream
   liquidOut.p = liquidIn.p;

   liquidOut.Ndot = -N_flow_delayed;
   liquidOut.h_outflow = inStream(liquidIn.h_outflow);
   liquidOut.x_outflow = inStream(liquidIn.x_outflow);
   Ndot_l=liquidIn.Ndot;

annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
            {100,100}}), graphics={
        Text(
          extent={{-100,-60},{100,-100}},
          lineColor={0,0,0},
          textString="%name"),
        Rectangle(extent={{-88,94},{92,-96}}, lineColor={0,0,255}),
        Line(
          points={{-66,-52},{-64,-18},{-58,10},{-38,38},{-12,54},{24,62},{70,62}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-66,-52},{-66,90}},
          color={0,0,0},
          smooth=Smooth.None,
          thickness=0.5,
          arrow={Arrow.None,Arrow.Filled}),
        Line(
          points={{-66,-52},{90,-52}},
          color={0,0,0},
          smooth=Smooth.None,
          thickness=0.5,
          arrow={Arrow.None,Arrow.Filled})}));
end Delay;
