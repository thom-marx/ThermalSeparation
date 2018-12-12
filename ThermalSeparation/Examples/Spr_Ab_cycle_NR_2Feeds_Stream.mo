within ThermalSeparation.Examples;
model Spr_Ab_cycle_NR_2Feeds_Stream
    //Lauftest: 1.2.
//Lauftest: 11.11.2011, 37s (Knudsen)
//Lauftest: 18.10.2011, 25 s (Laptop, 1. Commit, 2. Commit), check luft nicht
//Lauftest: 14.10.2011, 21 s (Kevin, 2. Commit), check luft nicht
//Lauftest: 13.10.2011, 24 s (altes BaseColumn), 25s (BaseColumn 4, Filmmodell: Equilibrium, 3. Commit)

ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid(
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O,
    x={0,0,0,0,0,0,1},
    T=323.15,
    Flow=0.03)            annotation (Placement(transformation(extent={{-18,70},
            {2,90}}, rotation=0)));

ThermalSeparation.Components.SourcesSinks.SinkLiquid sinkLiquid(
                                                     redeclare package Medium =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O)
                      annotation (Placement(transformation(extent={{-10,-10},{
            10,10}},
                   rotation=270,
        origin={34,-62})));

ThermalSeparation.Components.SourcesSinks.SinkGas sinkGas(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_SO2_HCl_HF_N2,
    p=149910)   annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-28,50})));

ThermalSeparation.Components.Columns.SprayColumn baseStage_GF1(
    mapping={{1,7},{2,2},{3,3},{4,4},{5,5},{6,6},{7,1}},
    inertVapour={false,false,false,false,false,false,false},
    redeclare package MediumVapour =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_SO2_HCl_HF_N2,
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O,
    nS=7,
    hasVapourFeed=false,
    n_elements=8,
    redeclare model Reaction = ThermalSeparation.Reaction.NoReaction,
    x_v_start_const={0.07,0.06,0.1,0.01,0.005,0.005,0.75},
    x_l_start_const={1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1 - 6e-5},
    hasLiquidFeed=true,
    stageLiquidFeed={6},
    redeclare model HeatTransferWall = ThermalSeparation.Wall.Adiabatic,
    redeclare record Geometry =
        ThermalSeparation.Geometry.StructuredPackedColumn.Geometry (
                                                          d=5),
    p_v_start_inlet=156000,
    p_v_start_outlet=152000)                 annotation (Placement(transformation(extent={{-34,-14},
            {14,32}}, rotation=0)));

  ThermalSeparation.Components.SourcesSinks.SourceGas
                                    sourceGas_Vdot(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_SO2_HCl_HF_N2,
    x={0.065,0.063,0.118,0.001,0.001,0.001,0.751},
    T=417.15)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
        origin={-28,-44})));

  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    offset=120,
    startTime=100,
    height=0) annotation (Placement(transformation(extent={{-68,-64},{-48,-44}},
          rotation=0)));
  ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink
    annotation (Placement(transformation(extent={{24,-10},{44,10}},
                                                                  rotation=0)));

  inner SystemTS systemTS
    annotation (Placement(transformation(extent={{-64,2},{-44,22}})));
  ThermalSeparation.Components.SourcesSinks.SplitLiquid_1p splitLiquid_x(
    redeclare package Medium =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O)
                annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={38,-28})));

  ThermalSeparation.Components.SourcesSinks.SplitLiquid_1p splitLiquid_x1(
    redeclare package Medium =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O)
                annotation (Placement(transformation(
        origin={58,34},
        extent={{-10,-10},{10,10}},
        rotation=0)));

  ThermalSeparation.Components.SourcesSinks.CombLiquid_x
                                       combLiquid_x(
    redeclare package Medium =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O,
    x_start={0,0,0,0,0,0,1})
    annotation (Placement(transformation(
        origin={84,56},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  Modelica.Blocks.Sources.Constant const(k=0.8)
    annotation (Placement(transformation(extent={{0,-40},{14,-26}})));
  Modelica.Blocks.Sources.Constant const1(k=0.727) annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={34,60})));
equation
  connect(ramp.y, sourceGas_Vdot.Flow_in) annotation (Line(
      points={{-47,-54},{-32,-54}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const.y, splitLiquid_x.split) annotation (Line(
      points={{14.7,-33},{26.35,-33},{26.35,-28},{32,-28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(splitLiquid_x1.split, const1.y) annotation (Line(
      points={{58,40},{58,60},{45,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(combLiquid_x.liquidPortIn2, sourceLiquid.liquidPortOut) annotation (
      Line(
      points={{94,61},{94,80},{3.4,80}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(baseStage_GF1.downStreamIn, splitLiquid_x1.liquidPortOut2)
    annotation (Line(
      points={{6.8,29.7},{26.4,29.7},{26.4,39},{48,39}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(baseStage_GF1.feedLiquid[1], splitLiquid_x1.liquidPortOut1)
    annotation (Line(
      points={{-30.16,8.08},{-30.16,14},{48,14},{48,29}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(baseStage_GF1.heatPort, ambientHeatSink.heatPort1) annotation (Line(
      points={{9.68,9},{18,9},{18,0},{22.4,0}},
      color={188,51,69},
      smooth=Smooth.None));
  connect(splitLiquid_x.liquidPortIn, baseStage_GF1.downStreamOut) annotation (
      Line(
      points={{38,-18},{6.8,-18},{6.8,-11.7}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(splitLiquid_x.liquidPortOut2, sinkLiquid.liquidPortIn) annotation (
      Line(
      points={{33,-38},{34,-38},{34,-52.4}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(splitLiquid_x.liquidPortOut1, combLiquid_x.liquidPortIn1) annotation (
     Line(
      points={{43,-38},{94,-38},{94,51}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(splitLiquid_x1.liquidPortIn, combLiquid_x.liquidPortOut) annotation (
      Line(
      points={{68,34},{68,56},{74,56}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(sourceGas_Vdot.gasPortOut, baseStage_GF1.upStreamIn) annotation (Line(
      points={{-28,-32.6},{-28,-11.7},{-26.8,-11.7}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(baseStage_GF1.upStreamOut, sinkGas.gasPortIn) annotation (Line(
      points={{-26.8,29.7},{-26.8,34.85},{-28,34.85},{-28,40.4}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                    graphics),
    experiment(StopTime=2000, Algorithm="Dassl"),
    experimentSetupOutput);
end Spr_Ab_cycle_NR_2Feeds_Stream;
