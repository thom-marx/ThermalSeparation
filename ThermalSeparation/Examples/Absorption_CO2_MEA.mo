within ThermalSeparation.Examples;
model Absorption_CO2_MEA "absorption of CO2 in MEA solution"
  //Lauftest: 1.2.
  //Lauftest: 25.1.
  //Lauftest: 12.1. (3s, Kevin)

ThermalSeparation.Components.SourcesSinks.SourceLiquid
                                              sourceLiquid(
    redeclare package MediumLiquid = ThermalSeparation.Media.H2O_CO2_MEA_Liq,
    x={0.868,0.022,0.11},
    Flow=1.1,
    T=313.15)             annotation (Placement(transformation(extent={{-10,-10},
            {10,10}}, rotation=270,
        origin={6,66})));
ThermalSeparation.Components.SourcesSinks.SinkLiquid
                                            sinkLiquid(redeclare package Medium =
        ThermalSeparation.Media.H2O_CO2_MEA_Liq)
                      annotation (Placement(transformation(extent={{-10,-10},{10,
            10}},  rotation=270,
        origin={6,-42})));

ThermalSeparation.Components.SourcesSinks.SinkGas
                                         sinkGas(redeclare package Medium =
        ThermalSeparation.Media.H2O_O2_CO2_N2_Vap, p=100000)
                annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-30,64})));

ThermalSeparation.Components.Columns.StructuredPackedColumn column(
    hasLiquidFeed=false,
    hasVapourFeed=false,
    mapping={{1,1},{3,2}},
    inertVapour={false,true,false,true},
    inertLiquid={false,false,true},
    redeclare package MediumVapour = ThermalSeparation.Media.H2O_O2_CO2_N2_Vap,
    redeclare package MediumLiquid = ThermalSeparation.Media.H2O_CO2_MEA_Liq,
    nS=2,
    redeclare model Reaction =
        ThermalSeparation.Reaction.ReactionEquilibrium.CO2_MEA,
    redeclare model InitOption =
        ThermalSeparation.Components.Columns.BaseClasses.Initialization.Init_T_xv_p_Ndot0,
    wettedInitial=false,
    eps_liq_start=0.05,
    redeclare model HeatTransferWall = ThermalSeparation.Wall.Adiabatic,
    redeclare record Geometry =
        ThermalSeparation.Geometry.StructuredPackedColumn.Mellapak250Y (d=15, H=
           15),
    T_l_profile=false,
    redeclare model PressureLoss =
        ThermalSeparation.PressureLoss.StructuredPackedColumn.Stichlmair,
    redeclare model BalanceEquations =
        ThermalSeparation.BalanceEquations.StructuredPackedColumn.NonEquilibrium.TwoPhaseVarState
        (redeclare model FilmModel =
            ThermalSeparation.FilmModel.StructuredPackedColumn.MS),
    redeclare model HomotopyMethod =
        ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.NoHomotopy,
    x_v_start_const={0.1841,0.0471,0.04,0.7768},
    x_l_start_const={0.85,0.05,0.1},
    n_elements=2,
    p_v_start_inlet=101000,
    p_v_start_outlet=100000,
    T_vapour_start=333.15,
    T_liquid_start=333.15,
    p_v_in(start=101000),
    redeclare model ThermoEquilibrium = PhaseEquilibrium.H2O_CO2_MEA)
                         annotation (Placement(transformation(extent={{-36,-6},
            {12,40}}, rotation=0)));

  ThermalSeparation.Components.SourcesSinks.SourceGas
                                             sourceGas_Vdot(
    redeclare package Medium = ThermalSeparation.Media.H2O_O2_CO2_N2_Vap,
    x={0.075,0.034,0.145,0.746},
    flowOption=ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowMdot,
    T=308.15)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
        origin={-28,-42})));

  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    startTime=100,
    height=0,
    offset=330)
              annotation (Placement(transformation(extent={{-68,-64},{-48,-44}},
          rotation=0)));
  ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink
    annotation (Placement(transformation(extent={{24,12},{36,24}},rotation=0)));
  inner SystemTS systemTS
    annotation (Placement(transformation(extent={{-100,80},{-80,100}})));

equation
  connect(ramp.y, sourceGas_Vdot.Flow_in) annotation (Line(
      points={{-47,-54},{-44,-54},{-44,-52},{-32,-52}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ambientHeatSink.heatPort1, column.heatPort) annotation (Line(
      points={{23.04,18},{18,18},{18,17},{8.16,17}},
      color={188,51,69},
      smooth=Smooth.None));
  connect(column.downStreamIn, sourceLiquid.liquidPortOut) annotation (Line(
      points={{4.8,37.7},{4.8,45.85},{6,45.85},{6,54.6}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(column.upStreamOut, sinkGas.gasPortIn) annotation (Line(
      points={{-28.8,37.7},{-28.8,45.85},{-30,45.85},{-30,54.4}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column.upStreamIn, sourceGas_Vdot.gasPortOut) annotation (Line(
      points={{-28.8,-3.7},{-28.8,-16.85},{-28,-16.85},{-28,-30.6}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column.downStreamOut, sinkLiquid.liquidPortIn) annotation (Line(
      points={{4.8,-3.7},{4.8,-18.85},{6,-18.85},{6,-32.4}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),
    experiment(StopTime=2000, Algorithm="Dassl"),
    experimentSetupOutput,
    Documentation(info="<html>
<p><ul>
<li> packed column</li>
<li> 3 components in vapour phase</li>
<li> 3 components in liquid phase</li>
<li> absorption</li>
</ul></p>
</html>"));
end Absorption_CO2_MEA;
