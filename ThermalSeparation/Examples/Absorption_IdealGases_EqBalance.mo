within ThermalSeparation.Examples;
model Absorption_IdealGases_EqBalance
  "absorption of ideal gases in water, using a packed column"
  //Lauftest: 1.2.
  //Lauftest: 25.1.
  //Lauftest: 12.1. (3s, Kevin)

ThermalSeparation.Components.SourcesSinks.SourceLiquid
                                              sourceLiquid(
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.CO2_water,
    T=323.15,
    x={0,1})              annotation (Placement(transformation(extent={{-10,-10},
            {10,10}}, rotation=270,
        origin={6,88})));
ThermalSeparation.Components.SourcesSinks.SinkLiquid
                                            sinkLiquid(redeclare package Medium =
        ThermalSeparation.Media.WaterBasedLiquid.CO2_water)
                      annotation (Placement(transformation(extent={{-10,-10},{10,
            10}},  rotation=270,
        origin={8,-30})));

ThermalSeparation.Components.SourcesSinks.SinkGas
                                         sinkGas(redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.CO2_H2O, p=149910)
                annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-30,86})));

  ThermalSeparation.Components.SourcesSinks.SourceGas
                                             sourceGas_Vdot(
    flowOption=ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowNdot,
    redeclare package Medium = ThermalSeparation.Media.IdealGasMixtures.CO2_H2O,
    x={0.85,0.15},
    T=417.15)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
        origin={-26,-28})));

  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    startTime=100,
    offset=450,
    height=0) annotation (Placement(transformation(extent={{-68,-64},{-48,-44}},
          rotation=0)));
  inner ThermalSeparation.SystemTS systemTS
    annotation (Placement(transformation(extent={{-74,-8},{-54,12}})));

ThermalSeparation.Components.Columns.StructuredPackedColumn column1(
    hasLiquidFeed=false,
    hasVapourFeed=false,
    redeclare model Reaction = ThermalSeparation.Reaction.NoReaction,
    redeclare record Geometry =
        ThermalSeparation.Geometry.StructuredPackedColumn.Geometry (H=10, d=4.9),
    x_l_start_const={1e-5,1 - 2e-5},
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.CO2_water,
    nS=2,
    x_v_start_const={0.96,0.04},
    mapping={{1,1},{2,2}},
    redeclare package MediumVapour =
        ThermalSeparation.Media.IdealGasMixtures.CO2_H2O,
    wettedInitial=false,
    eps_liq_start=0.3,
    redeclare model ThermoEquilibrium =
        ThermalSeparation.PhaseEquilibrium.IdealGasIdealLiquid,
    redeclare model InitOption =
        ThermalSeparation.Components.Columns.BaseClasses.Initialization.Init_EquilibriumBalancing,
    n_elements=3,
    redeclare model BalanceEquations =
        ThermalSeparation.BalanceEquations.StructuredPackedColumn.Equilibrium.TwoPhaseVarState
        (redeclare model FilmModel =
            ThermalSeparation.FilmModel.StructuredPackedColumn.EffectiveDiff),
    p_v_start_inlet=156000,
    p_v_start_outlet=152000)
                         annotation (Placement(transformation(extent={{-34,18},
            {14,64}}, rotation=0)));

  ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink1
    annotation (Placement(transformation(extent={{36,30},{56,50}},rotation=0)));
equation
  connect(ramp.y, sourceGas_Vdot.Flow_in) annotation (Line(
      points={{-47,-54},{-44,-54},{-44,-38},{-30,-38}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(column1.downStreamIn, sourceLiquid.liquidPortOut) annotation (Line(
      points={{6.8,61.7},{6.8,68.85},{6,68.85},{6,76.6}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(column1.upStreamOut, sinkGas.gasPortIn) annotation (Line(
      points={{-26.8,61.7},{-26.8,67.85},{-30,67.85},{-30,76.4}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column1.upStreamIn, sourceGas_Vdot.gasPortOut) annotation (Line(
      points={{-26.8,20.3},{-26.8,1.15},{-26,1.15},{-26,-16.6}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column1.downStreamOut, sinkLiquid.liquidPortIn) annotation (Line(
      points={{6.8,20.3},{6.8,0.15},{8,0.15},{8,-20.4}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(column1.heatPort, ambientHeatSink1.heatPort1) annotation (Line(
      points={{10.16,41},{22.08,41},{22.08,40},{34.4,40}},
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
end Absorption_IdealGases_EqBalance;
