within ThermalSeparation.Examples;
model Absorption_IdealGases_Random
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
        origin={10,62})));
ThermalSeparation.Components.SourcesSinks.SinkLiquid
                                            sinkLiquid(         redeclare package Medium =
                       ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O)
                      annotation (Placement(transformation(extent={{-10,-10},{10,
            10}},  rotation=270,
        origin={10,-54})));

ThermalSeparation.Components.SourcesSinks.SinkGas
                                         sinkGas(         redeclare package Medium =
               ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2, p=149910)
                annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-28,58})));

  ThermalSeparation.Components.SourcesSinks.SourceGas
                                             sourceGas_Vdot(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2,
    x={0.85,0.05,0.1},
    flowOption=ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowNdot,
    T=417.15)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
        origin={-24,-54})));

  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    startTime=100,
    height=0,
    offset=450)
              annotation (Placement(transformation(extent={{-66,-90},{-46,-70}},
          rotation=0)));
  inner SystemTS systemTS
    annotation (Placement(transformation(extent={{-100,80},{-80,100}})));

Components.Columns.RandomPackedColumn                       column1(
    hasLiquidFeed=false,
    hasVapourFeed=false,
    redeclare package MediumVapour =
        ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2,
    nS=3,
    mapping={{1,1},{2,3},{3,2}},
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O,
    redeclare model Reaction = ThermalSeparation.Reaction.NoReaction,
    eps_liq_start=0.02,
    x_l_start_const={1e-5,1e-5,1 - 2e-5},
    x_v_start_const={0.96,0.02,0.02},
    redeclare record Geometry =
        ThermalSeparation.Geometry.StructuredPackedColumn.Geometry (H=10, d=4.9),
    redeclare model ThermoEquilibrium =
        ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid (factor_K=
           {0.8,0.8,0.8}),
    wettedInitial=true,
    redeclare model PressureLoss =
        ThermalSeparation.PressureLoss.RandomPackedColumn.Particlemodel,
    T_l_profile=false,
    redeclare model HomotopyMethod =
        ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.NoHomotopy,
    n_elements=5,
    p_v_start_inlet=156000,
    p_v_start_outlet=152000,
    T_vapour_start=328.15,
    T_liquid_start=328.15,
    redeclare model BalanceEquations =
        ThermalSeparation.BalanceEquations.RandomPackedColumn.NonEquilibrium.TwoPhaseVarState (
         redeclare model FilmModel =
            ThermalSeparation.FilmModel.RandomPackedColumn.MS (redeclare model StateSelection =
                ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.StateSelection1)))
                         annotation (Placement(transformation(extent={{-32,-8},
            {16,38}}, rotation=0)));

  ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink1
    annotation (Placement(transformation(extent={{26,5},{46,25}}, rotation=0)));
equation
  connect(ramp.y, sourceGas_Vdot.Flow_in) annotation (Line(
      points={{-45,-80},{-42,-80},{-42,-64},{-28,-64}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(column1.upStreamOut, sinkGas.gasPortIn) annotation (Line(
      points={{-24.8,35.7},{-24.8,41.85},{-28,41.85},{-28,48.4}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column1.downStreamIn, sourceLiquid.liquidPortOut) annotation (Line(
      points={{8.8,35.7},{8.8,42.85},{10,42.85},{10,50.6}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(column1.heatPort, ambientHeatSink1.heatPort1) annotation (Line(
      points={{12.16,15},{24.4,15}},
      color={188,51,69},
      smooth=Smooth.None));
  connect(sinkLiquid.liquidPortIn, column1.downStreamOut) annotation (Line(
      points={{10,-44.4},{10,-5.7},{8.8,-5.7}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(sourceGas_Vdot.gasPortOut, column1.upStreamIn) annotation (Line(
      points={{-24,-42.6},{-24,-5.7},{-24.8,-5.7}},
      color={255,127,39},
      thickness=1,
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
end Absorption_IdealGases_Random;
