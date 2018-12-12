within ThermalSeparation.Components.Sensors;
model TemperatureSensor
replaceable package Medium = Media.BaseMediumLiquid annotation(choicesAllMatching);
final parameter Integer nSL = Medium.nSubstance;

Medium.BaseProperties medium(p=liquidIn.p, x=inStream(liquidIn.x_outflow),T=T_sensor,h=h_l);
SI.Temperature T_sensor;
SI.MolarEnthalpy h_l;

   ThermalSeparation.Interfaces.LiquidPortIn
                           liquidIn(redeclare package Medium=Medium)
                                        annotation (Placement(transformation(
          extent={{-110,-10},{-90,10}},rotation=0), iconTransformation(extent={{-110,
            -10},{-90,10}})));

   ThermalSeparation.Interfaces.LiquidPortOut
                            liquidOut(redeclare package Medium=Medium)
                                          annotation (Placement(transformation(
          extent={{90,-10},{110,10}},
                                   rotation=0), iconTransformation(extent={{90,-10},
            {110,10}})));

  Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,102}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-2,100})));

equation
  medium.h = inStream(liquidIn.h_outflow);
   liquidIn.p = liquidOut.p;

   liquidIn.Ndot = -liquidOut.Ndot;
   inStream(liquidIn.h_outflow) = liquidOut.h_outflow;
   inStream(liquidIn.x_outflow) = liquidOut.x_outflow;

      inStream(liquidOut.h_outflow) = liquidIn.h_outflow;
   inStream(liquidOut.x_outflow) = liquidIn.x_outflow;

   y = medium.T;
  annotation (Diagram(graphics), Icon(graphics={
        Polygon(
          points={{18,-70},{58,-85},{18,-100},{18,-70}},
          lineColor={0,128,255},
          smooth=Smooth.None,
          fillColor={0,128,255},
          fillPattern=FillPattern.Solid,
          visible=showDesignFlowDirection),
        Polygon(
          points={{18,-75},{48,-85},{18,-95},{18,-75}},
          lineColor={255,255,255},
          smooth=Smooth.None,
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          visible=allowFlowReversal),
        Line(
          points={{53,-85},{-62,-85}},
          color={0,128,255},
          smooth=Smooth.None,
          visible=showDesignFlowDirection),
        Text(
          extent={{-139,-104},{161,-144}},
          lineColor={0,0,255},
          textString="%name"),
        Line(points={{-2,100},{-2,50}},
                                      color={0,0,127}),
        Line(points={{-94,0},{98,0}},  color={0,128,255}),
        Ellipse(
          extent={{-22,-68},{18,-30}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-14,50},{10,-34}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-14,50},{-14,70},{-12,76},{-8,78},{-2,80},{4,78},{8,76},{10,70},
              {10,50},{-14,50}},
          lineColor={0,0,0},
          lineThickness=0.5),
        Line(
          points={{-14,50},{-14,-35}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{10,50},{10,-34}},
          color={0,0,0},
          thickness=0.5),
        Line(points={{-42,-10},{-14,-10}}, color={0,0,0}),
        Line(points={{-42,20},{-14,20}}, color={0,0,0}),
        Line(points={{-42,50},{-14,50}}, color={0,0,0}),
        Text(
          extent={{92,122},{-2,92}},
          lineColor={0,0,0},
          textString="T")}));
end TemperatureSensor;
