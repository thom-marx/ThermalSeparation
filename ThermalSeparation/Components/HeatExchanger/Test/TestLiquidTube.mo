within ThermalSeparation.Components.HeatExchanger.Test;
model TestLiquidTube

  LiquidTube liquidTube(
    height_in=0,
    zeta=0.02,
    redeclare package MediumLiquid = 
        ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS_WaterProp,
    c_l_start={0,0,0}) 
    annotation (Placement(transformation(extent={{-38,-40},{14,14}})));
  ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid_x(
    redeclare package MediumLiquid = 
        ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS_WaterProp,
    x={0,1,0},
    Flow=1e-3) annotation (Placement(transformation(extent={{-36,42},{0,78}})));
  SourcesSinks.SinkLiquid sinkLiquid(           redeclare package Medium = 
        ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS_DesorberGamma)
    annotation (Placement(transformation(extent={{30,-98},{64,-66}})));
  inner SystemTS systemTS 
    annotation (Placement(transformation(extent={{-72,28},{-52,48}})));
equation
  connect(sourceLiquid_x.liquidPortOut, liquidTube.liquidIn) annotation (Line(
      points={{-18,44.52},{-18,27.72},{-12,27.72},{-12,13.46}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  connect(sinkLiquid.liquidPortIn, liquidTube.liquidOut) annotation (Line(
      points={{47,-70.48},{16.5,-70.48},{16.5,-39.46},{-12,-39.46}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  annotation (Diagram(graphics));
end TestLiquidTube;
