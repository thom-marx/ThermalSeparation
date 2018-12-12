within ThermalSeparation.Components.Pumps.Test;
model idealPumpTest
//Lauftest: 1.7.2011
  idealPump idealPump1
    annotation (Placement(transformation(extent={{-20,20},{0,40}})));
  SourcesSinks.SinkLiquid sinkLiquid_p1(
    redeclare package Medium = ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
    p=600000)
    annotation (Placement(transformation(extent={{20,-8},{40,12}})));

  ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid_x_p(
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
    x={0,1},
    use_Flow=true,
    fixed_pressure=true,
    p_fixed=300000) annotation (Placement(transformation(extent={{-62,30},{-42,50}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=0.0003)
    annotation (Placement(transformation(extent={{-92,52},{-72,72}})));
  inner SystemTS systemTS
    annotation (Placement(transformation(extent={{-92,8},{-72,28}})));
equation
  connect(sinkLiquid_p1.liquidPortIn, idealPump1.liquidOut) annotation (Line(
      points={{30,9.2},{18,9.2},{18,36},{-2.4,36}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  connect(sourceLiquid_x_p.liquidPortOut, idealPump1.liquidIn) annotation (Line(
      points={{-52,31.4},{-35.1,31.4},{-35.1,30},{-17.4,30}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  connect(realExpression.y, sourceLiquid_x_p.Flow_In) annotation (Line(
      points={{-71,62},{-71,53},{-60.4,53},{-60.4,42.8}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
            {100,100}}), graphics));
end idealPumpTest;
