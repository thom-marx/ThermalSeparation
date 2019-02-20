within ThermalSeparation.Examples;
model ProPen_EQ_SteadyState
  "separation of alkanes, using a packed column, steady state balance equations"
  //Lauftest: 1.2.
  //Lauftest: 25.1.
  //Lauftest: 12.1. (3s, Kevin)

ThermalSeparation.Components.SourcesSinks.SourceLiquid
                                              sourceLiquid(
    redeclare package MediumLiquid =
        ThermalSeparation.Media.Propane_Pentane_Liq,
    x={0.3,0.7},
    Flow=1,
    T=298.15)             annotation (Placement(transformation(extent={{-10,-10},
            {10,10}}, rotation=270,
        origin={6,90})));
ThermalSeparation.Components.SourcesSinks.SinkLiquid
                                            sinkLiquid(redeclare package Medium =
        ThermalSeparation.Media.Propane_Pentane_Liq)
                      annotation (Placement(transformation(extent={{-10,-10},{
            10,10}},
                   rotation=270,
        origin={8,-32})));

ThermalSeparation.Components.SourcesSinks.SinkGas
                                         sinkGas(redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.Propane_Pentane, p=300000)
                annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-30,88})));

  ThermalSeparation.Components.SourcesSinks.SourceGas
                                             sourceGas_Vdot(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.Propane_Pentane,
    x={0.5,0.5},
    T=323.15,
    flowOption=ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowNdot)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
        origin={-28,-32})));

  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    startTime=100,
    height=0,
    offset=450)
              annotation (Placement(transformation(extent={{-68,-64},{-48,-44}},
          rotation=0)));
  inner ThermalSeparation.SystemTS systemTS
    annotation (Placement(transformation(extent={{-74,-8},{-54,12}})));

ThermalSeparation.Components.Columns.StructuredPackedColumn column1(
    hasLiquidFeed=false,
    hasVapourFeed=false,
    redeclare model Reaction = ThermalSeparation.Reaction.NoReaction,
    nS=2,
    mapping={{1,1},{2,2}},
    eps_liq_start=0.3,
    redeclare package MediumVapour =
        ThermalSeparation.Media.IdealGasMixtures.Propane_Pentane,
    redeclare package MediumLiquid =
        ThermalSeparation.Media.Propane_Pentane_Liq,
    redeclare record Geometry =
        ThermalSeparation.Geometry.StructuredPackedColumn.Geometry (H=10, d=4.9),
    x_l_start_const={0.5,0.5},
    x_v_start_const={0.3,0.7},
    wettedInitial=false,
    redeclare model InitOption =
        ThermalSeparation.Components.Columns.BaseClasses.Initialization.None,
    redeclare model ThermoEquilibrium =
        ThermalSeparation.PhaseEquilibrium.IdealGasIdealLiquid,
    redeclare model PressureLoss =
        ThermalSeparation.PressureLoss.StructuredPackedColumn.NominalLinear,
    redeclare model HeatTransferWall = ThermalSeparation.Wall.Adiabatic (stat=
            true),
    redeclare model BalanceEquations =
        ThermalSeparation.BalanceEquations.StructuredPackedColumn.Equilibrium.TwoPhaseSteadyState,
    n_elements=1,
    p_v_start_inlet=302000,
    p_v_start_outlet=300000,
    T_vapour_start=303.15,
    T_liquid_start=303.15)
                         annotation (Placement(transformation(extent={{-34,18},
            {14,64}}, rotation=0)));

  ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink1
    annotation (Placement(transformation(extent={{40,31},{60,51}},rotation=0)));
equation
  connect(ramp.y, sourceGas_Vdot.Flow_in) annotation (Line(
      points={{-47,-54},{-44,-54},{-44,-42},{-32,-42}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sinkGas.gasPortIn, column1.upStreamOut) annotation (Line(
      points={{-30,78.4},{-30,61.7},{-26.8,61.7}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column1.downStreamIn, sourceLiquid.liquidPortOut) annotation (Line(
      points={{6.8,61.7},{6.8,68.85},{6,68.85},{6,78.6}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(sinkLiquid.liquidPortIn, column1.downStreamOut) annotation (Line(
      points={{8,-22.4},{6.8,-22.4},{6.8,20.3}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(sourceGas_Vdot.gasPortOut, column1.upStreamIn) annotation (Line(
      points={{-28,-20.6},{-28,-4},{-26.8,-4},{-26.8,20.3}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column1.heatPort, ambientHeatSink1.heatPort1) annotation (Line(
      points={{10.16,41},{38.4,41}},
      color={188,51,69},
      smooth=Smooth.None));
annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),
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
end ProPen_EQ_SteadyState;
