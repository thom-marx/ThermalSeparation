within ThermalSeparation.Components.HeatExchanger.Test;
model counterFlowHeatExchangerH2OTest

  counterFlowHeatExchangerH2O counterFlowHeatExchanger1(
    plateMetalMass=20,
    volumeHotLiquid(displayUnit="l") = 0.001,
    volumeColdLiquid(displayUnit="l") = 0.001,
    mFlowColdLiquid_nom=0.42,
    surfaceArea=0.4,
    alphaHot=5400,
    alphaCold=5400,
    dpHotLiquid=50000,
    dpColdLiquid=50000,
    ThotLiquidIn_start=329.15,
    TcoldLiquidIn_start=298.15,
    ThotLiquid_start=313.15,
    TcoldLiquid_start=308.15) 
    annotation (Placement(transformation(extent={{-30,-8},{8,48}})));
  SourcesSinks.SinkLiquid sinkLiquid_p1(redeclare package Medium = 
        ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O, p=200000) 
    annotation (Placement(transformation(extent={{-70,-28},{-50,-8}})));
  ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid_x1(
    redeclare package MediumLiquid = 
        ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
    x={0.1,0.9},
    T=329.15,
    Flow=3.1e-4) 
    annotation (Placement(transformation(extent={{30,-10},{50,10}})));
  Modelica.Blocks.Sources.Ramp ramp(
    startTime=500,
    duration=50,
    height=0.1,
    offset=0.42) 
    annotation (Placement(transformation(extent={{-100,58},{-80,78}})));
  Modelica.Fluid.Sources.MassFlowSource_T boundary(
    use_m_flow_in=true,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    nPorts=1,
    T=298.15) annotation (Placement(transformation(extent={{-70,48},{-50,68}})));
  Modelica.Fluid.Sources.FixedBoundary boundary1(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    nPorts=1,
    p=200000,
    T=324.15) annotation (Placement(transformation(extent={{56,28},{36,48}})));
  inner SystemTS systemTS 
    annotation (Placement(transformation(extent={{-82,18},{-62,38}})));
equation
  connect(counterFlowHeatExchanger1.hotLiquidOut, sinkLiquid_p1.liquidPortIn) 
    annotation (Line(
      points={{-30,0.4},{-60,0.4},{-60,-10.8}},
      smooth=Smooth.None,
      color={0,0,0}));
  connect(sourceLiquid_x1.liquidPortOut, counterFlowHeatExchanger1.hotLiquidIn) 
    annotation (Line(
      points={{40,-8.6},{40,-20},{20,-20},{20,0.4},{8,0.4}},
      smooth=Smooth.None,
      color={0,0,0}));
  connect(boundary.ports[1], counterFlowHeatExchanger1.port_a) annotation (Line(
      points={{-50,58},{-40,58},{-40,37.36},{-30,37.36}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(counterFlowHeatExchanger1.port_b, boundary1.ports[1]) annotation (
      Line(
      points={{8.38,38.48},{22.19,38.48},{22.19,38},{36,38}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(ramp.y, boundary.m_flow_in) annotation (Line(
      points={{-79,68},{-74,68},{-74,66},{-70,66}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics));
end counterFlowHeatExchangerH2OTest;
