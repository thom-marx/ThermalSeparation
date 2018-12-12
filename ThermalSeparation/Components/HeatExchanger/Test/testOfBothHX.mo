within ThermalSeparation.Components.HeatExchanger.Test;
model testOfBothHX

  counterFlowHeatExchangerH2O W103(
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
    annotation (Placement(transformation(extent={{-24,-96},{14,-40}})));
  SourcesSinks.SinkLiquid sinkLiquid_p1(
    use_p=true,
    redeclare package Medium = ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
    p=200000) 
    annotation (Placement(transformation(extent={{-74,-66},{-52,-84}})));

  Modelica.Blocks.Sources.Ramp coolingWater(
    startTime=500,
    duration=50,
    height=0.1,
    offset=0.42) 
    annotation (Placement(transformation(extent={{-98,-12},{-62,24}})));
  Modelica.Fluid.Sources.MassFlowSource_T boundary(
    use_m_flow_in=true,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    nPorts=1,
    T=298.15) 
    annotation (Placement(transformation(extent={{-54,-38},{-34,-18}})));
  Modelica.Fluid.Sources.FixedBoundary boundary1(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    nPorts=1,
    p=200000,
    T=324.15) annotation (Placement(transformation(extent={{50,-62},{30,-42}})));
  counterFlowHeatExchanger W102(volumeHotLiquid(displayUnit="l"),
      volumeColdLiquid(displayUnit="l")) 
    annotation (Placement(transformation(extent={{110,-96},{148,-40}})));
  SourcesSinks.SinkLiquid sinkLiquid_p(
    redeclare package Medium = ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
    p=400000) 
    annotation (Placement(transformation(extent={{170,-54},{190,-34}})));

  ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid_x2(
    redeclare package MediumLiquid = 
        ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
    x={0.1,0.9},
    T=377.15,
    Flow=3.1e-4) 
    annotation (Placement(transformation(extent={{170,-94},{190,-74}})));
  ThermalSeparation.Components.SourcesSinks.SourceLiquid Asborber(
    redeclare package MediumLiquid = 
        ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
    x={0.1,0.9},
    T=324.15,
    Flow=3.05e-4) 
    annotation (Placement(transformation(extent={{38,4},{116,68}})));
  inner SystemTS systemTS 
    annotation (Placement(transformation(extent={{-48,40},{-28,60}})));
equation
  connect(W103.hotLiquidOut, sinkLiquid_p1.liquidPortIn) 
    annotation (Line(
      points={{-24,-87.6},{-63,-87.6},{-63,-81.48}},
      smooth=Smooth.None,
      color={0,0,0}));
  connect(boundary.ports[1], W103.port_a) annotation (Line(
      points={{-34,-28},{-34,-50.64},{-24,-50.64}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(W103.port_b, boundary1.ports[1]) annotation (Line(
      points={{14.38,-49.52},{28.19,-49.52},{28.19,-52},{30,-52}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(coolingWater.y, boundary.m_flow_in) annotation (Line(
      points={{-60.2,6},{-68,6},{-68,-20},{-54,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(W102.coldLiquidOut, sinkLiquid_p.liquidPortIn) 
    annotation (Line(
      points={{142.68,-50.64},{142.68,-24},{180,-24},{180,-36.8}},
      smooth=Smooth.None,
      color={0,0,0}));
  connect(Asborber.liquidPortOut, W102.coldLiquidIn) 
    annotation (Line(
      points={{77,8.48},{77,-50.08},{115.32,-50.08}},
      smooth=Smooth.None,
      color={0,0,0}));
  connect(sourceLiquid_x2.liquidPortOut, W102.hotLiquidIn) annotation (Line(
      points={{180,-92.6},{163.9,-92.6},{163.9,-85.36},{142.3,-85.36}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  connect(W102.hotLiquidOut, W103.hotLiquidIn) annotation (Line(
      points={{115.32,-85.36},{64.66,-85.36},{64.66,-87.6},{14,-87.6}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{200,100}}), graphics), Icon(coordinateSystem(
          preserveAspectRatio=true, extent={{-100,-100},{200,100}})));
end testOfBothHX;
