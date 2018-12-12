within ThermalSeparation.Components.HeatExchanger.Test;
model counterFlowHeatExchangerTest

  counterFlowHeatExchanger counterFlowHeatExchanger1(volumeHotLiquid(
        displayUnit="l"), volumeColdLiquid(displayUnit="l")) 
    annotation (Placement(transformation(extent={{-30,-8},{8,48}})));
  SourcesSinks.SinkLiquid sinkLiquid_p(
    use_p=true,
    redeclare package Medium = ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
    p=400000) 
    annotation (Placement(transformation(extent={{30,30},{50,50}})));

  SourcesSinks.SinkLiquid sinkLiquid_p1(
    redeclare package Medium = ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
    p=200000) 
    annotation (Placement(transformation(extent={{-70,-28},{-50,-8}})));

  ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid_x1(
    redeclare package MediumLiquid = 
        ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
    x={0.1,0.9},
    Flow=3.1e-4,
    T=377.15) 
    annotation (Placement(transformation(extent={{30,-10},{50,10}})));
  Modelica.Blocks.Sources.Ramp ramp(
    startTime=500,
    height=1e-4,
    duration=50,
    offset=3.046e-4) 
    annotation (Placement(transformation(extent={{-100,50},{-80,70}})));
  ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid_x2(
    redeclare package MediumLiquid = 
        ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
    x={0.1,0.9},
    T=324.15,
    use_Flow=true) 
    annotation (Placement(transformation(extent={{-70,46},{-50,66}})));
  inner SystemTS systemTS 
    annotation (Placement(transformation(extent={{-82,18},{-62,38}})));
equation
  connect(counterFlowHeatExchanger1.coldLiquidOut, sinkLiquid_p.liquidPortIn) 
    annotation (Line(
      points={{2.68,37.36},{2.68,60},{40,60},{40,47.2}},
      smooth=Smooth.None,
      color={0,0,0}));
  connect(counterFlowHeatExchanger1.hotLiquidOut, sinkLiquid_p1.liquidPortIn) 
    annotation (Line(
      points={{-24.68,2.64},{-60,2.64},{-60,-10.8}},
      smooth=Smooth.None,
      color={0,0,0}));
  connect(sourceLiquid_x1.liquidPortOut, counterFlowHeatExchanger1.hotLiquidIn) 
    annotation (Line(
      points={{40,-8.6},{40,-20},{20,-20},{20,2.64},{2.3,2.64}},
      smooth=Smooth.None,
      color={0,0,0}));
  connect(sourceLiquid_x2.liquidPortOut, counterFlowHeatExchanger1.coldLiquidIn) 
    annotation (Line(
      points={{-60,47.4},{-60,37.92},{-24.68,37.92}},
      smooth=Smooth.None,
      color={0,0,0}));
  connect(ramp.y, sourceLiquid_x2.Flow_In) annotation (Line(
      points={{-79,60},{-73.7,60},{-73.7,58.8},{-68.4,58.8}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics), Commands(file=
          "Components/HeatExchanger/Test/plotCounterFlowHXtest.mos"
        "plotCounterFlowHXtest"));
end counterFlowHeatExchangerTest;
