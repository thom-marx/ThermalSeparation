within ThermalSeparation.Examples;
model Absorption_IdealGases_4components
  "absorption of ideal gases (4 components), packed absorber"
      //Lauftest: 1.2.
  //Lauftest: 20.11.2011, 15s (Laptop)

ThermalSeparation.Components.SourcesSinks.SourceLiquid
                                              sourceLiquid(
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_H2O,
    x={0,0,0,1},
    T=323.15,
    Flow=0.04)             annotation (Placement(transformation(extent={{-10,-10},
            {10,10}}, rotation=270,
        origin={6,54})));
ThermalSeparation.Components.SourcesSinks.SinkLiquid
                                            sinkLiquid(         redeclare
      package Medium = ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_H2O)
                      annotation (Placement(transformation(extent={{-10,-10},{
            10,10}},
                   rotation=270,
        origin={8,-42})));

ThermalSeparation.Components.SourcesSinks.SinkGas
                                         sinkGas(         redeclare package
      Medium = ThermalSeparation.Media.IdealGasMixtures.N2_H2O_CO2_O2, p=149910)
                annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-28,52})));

ThermalSeparation.Components.Columns.StructuredPackedColumn column(
    hasLiquidFeed=false,
    hasVapourFeed=false,
    redeclare model Reaction = ThermalSeparation.Reaction.NoReaction,
    eps_liq_start=0.02,
    wettedInitial=true,
    redeclare record Geometry =
        ThermalSeparation.Geometry.StructuredPackedColumn.Geometry (H=10, d=4.9),
    x_l_start_const={1e-5,1e-5,1e-5,1 - 3e-5},
    mapping={{1,4},{2,2},{3,3},{4,1}},
    redeclare package MediumVapour =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_N2,
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_H2O,
    nS=4,
    x_v_start_const={0.001,0.001,0.001,0.997},
    p_v_start_inlet=156000,
    p_v_start_outlet=152000,
    n_elements=5)        annotation (Placement(transformation(extent={{-34,-14},
            {14,32}}, rotation=0)));

//    p_v_start={1.6e5,1.55e5},
  ThermalSeparation.Components.SourcesSinks.SourceGas
                                             sourceGas_Vdot(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.N2_H2O_CO2_O2,
    T=417.15,
    x={0.85,0.05,0.05,0.05})
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
        origin={-26,-42})));

  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    startTime=100,
    height=0,
    offset=10)
              annotation (Placement(transformation(extent={{-68,-64},{-48,-44}},
          rotation=0)));
  ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink
    annotation (Placement(transformation(extent={{24,-2},{44,18}},rotation=0)));
  inner SystemTS systemTS
    annotation (Placement(transformation(extent={{-74,-8},{-54,12}})));

equation
  connect(ramp.y, sourceGas_Vdot.Flow_in) annotation (Line(
      points={{-47,-54},{-44,-54},{-44,-52},{-30,-52}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ambientHeatSink.heatPort1, column.heatPort) annotation (Line(
      points={{22.4,8},{16,8},{16,9},{10.16,9}},
      color={188,51,69},
      smooth=Smooth.None));
  connect(column.downStreamIn, sourceLiquid.liquidPortOut) annotation (Line(
      points={{6.8,29.7},{6.8,34.85},{6,34.85},{6,42.6}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(sinkGas.gasPortIn, column.upStreamOut) annotation (Line(
      points={{-28,42.4},{-28,29.7},{-26.8,29.7}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column.upStreamIn, sourceGas_Vdot.gasPortOut) annotation (Line(
      points={{-26.8,-11.7},{-26.8,-20.85},{-26,-20.85},{-26,-30.6}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column.downStreamOut, sinkLiquid.liquidPortIn) annotation (Line(
      points={{6.8,-11.7},{6.8,-22.85},{8,-22.85},{8,-32.4}},
      color={153,217,234},
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
end Absorption_IdealGases_4components;
