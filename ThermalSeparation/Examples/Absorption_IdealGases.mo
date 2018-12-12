within ThermalSeparation.Examples;
model Absorption_IdealGases
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
        origin={6,92})));
ThermalSeparation.Components.SourcesSinks.SinkLiquid
                                            sinkLiquid(         redeclare
      package Medium = ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O)
                      annotation (Placement(transformation(extent={{-10,-10},{
            10,10}},
                   rotation=270,
        origin={7,-30})));

ThermalSeparation.Components.SourcesSinks.SinkGas
                                         sinkGas(         redeclare package
      Medium = ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2, p=149910)
                annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-26,90})));

  ThermalSeparation.Components.SourcesSinks.SourceGas
                                             sourceGas_Vdot(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2,
    x={0.85,0.05,0.1},
    flowOption=ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowNdot,
    T=417.15)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
        origin={-27,-28})));

  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    startTime=100,
    height=0,
    offset=450)
              annotation (Placement(transformation(extent={{-62,-41},{-56,-35}},
          rotation=0)));
  inner SystemTS systemTS
    annotation (Placement(transformation(extent={{-100,80},{-80,100}})));

ThermalSeparation.Components.Columns.StructuredPackedColumn column1(
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
    wettedInitial=true,
    x_l_start_const={1e-5,1e-5,1 - 2e-5},
    x_v_start_const={0.96,0.02,0.02},
    redeclare record Geometry =
        ThermalSeparation.Geometry.StructuredPackedColumn.Geometry (H=10, d=4.9),
    redeclare model ThermoEquilibrium =
        ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid (factor_K=
           {0.8,0.8,0.8}),
    n_elements=10,
    p_v_start_inlet=156000,
    p_v_start_outlet=152000,
    redeclare model HomotopyMethod =
        ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.NoHomotopy)
                         annotation (Placement(transformation(extent={{-32,18},{
            16,64}},  rotation=0)));

  ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink1
    annotation (Placement(transformation(extent={{36,33},{48,50}},rotation=0)));
equation
  connect(sourceGas_Vdot.gasPortOut, column1.upStreamIn) annotation (Line(
      points={{-27,-16.6},{-27,20.3},{-24.8,20.3}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column1.downStreamOut, sinkLiquid.liquidPortIn) annotation (Line(
      points={{8.8,20.3},{8.8,-0.85},{7,-0.85},{7,-20.4}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(sourceLiquid.liquidPortOut, column1.downStreamIn) annotation (Line(
      points={{6,80.6},{6,61.7},{8.8,61.7}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(sinkGas.gasPortIn, column1.upStreamOut) annotation (Line(
      points={{-26,80.4},{-26,61.7},{-24.8,61.7}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column1.heatPort, ambientHeatSink1.heatPort1) annotation (Line(
      points={{12.16,41},{24.08,41},{24.08,41.5},{35.04,41.5}},
      color={188,51,69},
      thickness=1,
      smooth=Smooth.None));
  connect(ramp.y, sourceGas_Vdot.Flow_in) annotation (Line(
      points={{-55.7,-38},{-31,-38}},
      color={0,0,127},
      smooth=Smooth.None));
annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                    graphics),
    Documentation(info="<html>
<p><ul>
<li> packed column</li>
<li> 3 components in vapour phase</li>
<li> 3 components in liquid phase</li>
<li> absorption</li>
</ul></p>
</html>"),
    experiment(StopTime=100));
end Absorption_IdealGases;
