within ThermalSeparation.Examples;
model Absorption_IdealGases_FeedStream2
  "absorption of ideal gases in water, using a packed column"
  //Lauftest: 1.2.
  //Lauftest: 25.1.
  //Lauftest: 12.1. (3s, Kevin)

ThermalSeparation.Components.SourcesSinks.SourceLiquid
                                              sourceLiquid(
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O,
    x={0,0,1},
    T=323.15)             annotation (Placement(transformation(extent={{-10,-10},
            {10,10}}, rotation=270,
        origin={8,58})));
ThermalSeparation.Components.SourcesSinks.SinkLiquid
                                            sinkLiquid(         redeclare package
              Medium = ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O)
                      annotation (Placement(transformation(extent={{-10,-10},{
            10,10}},
                   rotation=270,
        origin={8,-38})));

ThermalSeparation.Components.SourcesSinks.SinkGas
                                         sinkGas(         redeclare package
      Medium = ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2, p=149910)
                annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-24,58})));

ThermalSeparation.Components.Columns.StructuredPackedColumn column(
    redeclare package MediumVapour =
        ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2,
    nS=3,
    mapping={{1,1},{2,3},{3,2}},
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O,
    redeclare model Reaction = ThermalSeparation.Reaction.NoReaction,
    eps_liq_start=0.02,
    wettedInitial=true,
    x_l_start_const={1e-5,1e-5,1 - 2e-5},
    x_v_start_const={0.96,0.02,0.02},
    redeclare record Geometry =
        ThermalSeparation.Geometry.StructuredPackedColumn.Geometry (H=10, d=4.9),
    redeclare model ThermoEquilibrium =
        ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid (factor_K=
           {0.8,0.8,0.8}),
    hasVapourFeed=true,
    hasLiquidFeed=true,
    n_elements=5,
    stageLiquidFeed={3},
    stageVapourFeed={2},
    nonFeed_stages_l={1,2,4,5},
    nonFeed_stages_v={1,3,4,5},
    p_v_start_inlet=156000,
    p_v_start_outlet=152000)
                         annotation (Placement(transformation(extent={{-32,-12},
            {16,34}}, rotation=0)));

  ThermalSeparation.Components.SourcesSinks.SourceGas
                                             sourceGas_Vdot(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2,
    x={0.85,0.05,0.1},
    flowOption=ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowNdot,
    T=417.15)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
        origin={-24,-38})));

  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    startTime=100,
    height=0,
    offset=450)
              annotation (Placement(transformation(extent={{-68,-64},{-48,-44}},
          rotation=0)));
  ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink
    annotation (Placement(transformation(extent={{22,-2},{42,18}},rotation=0)));
  inner SystemTS systemTS
    annotation (Placement(transformation(extent={{-100,80},{-80,100}})));

  ThermalSeparation.Components.SourcesSinks.SourceGas
                                             sourceGas_Vdot1(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2,
    x={0.85,0.05,0.1},
    flowOption=ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowNdot,
    use_Flow=false,
    T=417.15)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=0,
        origin={-48,30})));
ThermalSeparation.Components.SourcesSinks.SourceLiquid
                                              sourceLiquid1(
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O,
    x={0,0,1},
    T=323.15)             annotation (Placement(transformation(extent={{-10,-10},
            {10,10}}, rotation=0,
        origin={-48,-6})));
equation
  connect(ramp.y, sourceGas_Vdot.Flow_in) annotation (Line(
      points={{-47,-54},{-44,-54},{-44,-48},{-28,-48}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sourceLiquid1.liquidPortOut, column.feedLiquid[1]) annotation (Line(
      points={{-36.6,-6},{-36.6,10.08},{-28.16,10.08}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(column.feedVapour[1], sourceGas_Vdot1.gasPortOut) annotation (Line(
      points={{-28.16,15.14},{-36.6,15.14},{-36.6,30}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(sourceGas_Vdot.gasPortOut, column.upStreamIn) annotation (Line(
      points={{-24,-26.6},{-24,-9.7},{-24.8,-9.7}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(sinkLiquid.liquidPortIn, column.downStreamOut) annotation (Line(
      points={{8,-28.4},{8,-9.7},{8.8,-9.7}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(column.downStreamIn, sourceLiquid.liquidPortOut) annotation (Line(
      points={{8.8,31.7},{8.8,38.85},{8,38.85},{8,46.6}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(column.upStreamOut, sinkGas.gasPortIn) annotation (Line(
      points={{-24.8,31.7},{-24.8,39.85},{-24,39.85},{-24,48.4}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(ambientHeatSink.heatPort1, column.heatPort) annotation (Line(
      points={{20.4,8},{16,8},{16,11},{12.16,11}},
      color={188,51,69},
      smooth=Smooth.None));
annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                    graphics),
    experiment(StopTime=60, Algorithm="Dassl"),
    experimentSetupOutput,
    Documentation(info="<html>
<p><ul>
<li> packed column</li>
<li> 3 components in vapour phase</li>
<li> 3 components in liquid phase</li>
<li> absorption</li>
</ul></p>
</html>"));
end Absorption_IdealGases_FeedStream2;
