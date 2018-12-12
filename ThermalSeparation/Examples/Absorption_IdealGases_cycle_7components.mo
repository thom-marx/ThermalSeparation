within ThermalSeparation.Examples;
model Absorption_IdealGases_cycle_7components
  "absorption of ideal gases in water (7 components), cylce, spray absorber"
      //Lauftest: 1.2.

ThermalSeparation.Components.SourcesSinks.SourceLiquid
                                              sourceLiquid(
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O,
    x={0,0,0,0,0,0,1},
    flowOption=ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowNdot,
    T=323.15,
    Flow=3.5e3)                    annotation (Placement(transformation(extent=
            {{2,46},{22,66}}, rotation=0)));

ThermalSeparation.Components.SourcesSinks.SinkLiquid
                                            sinkLiquid(
                                                    redeclare package Medium =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O)
                      annotation (Placement(transformation(extent={{-10,-10},{
            10,10}},
                   rotation=270,
        origin={20,-62})));

ThermalSeparation.Components.SourcesSinks.SinkGas
                                         sinkGas(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_SO2_HCl_HF_N2,
    p=149910)   annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-24,52})));

ThermalSeparation.Components.Columns.SprayColumn baseStage_GF1(
    hasLiquidFeed=false,
    hasVapourFeed=false,
    redeclare model Reaction = ThermalSeparation.Reaction.NoReaction,
    redeclare model HeatTransferWall =
        ThermalSeparation.Wall.Adiabatic,
    nS=7,
    x_l_start_const={1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1 - 6e-5},
    inertVapour={false,false,false,false,false,false,false},
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O,
    x_v_start_const={0.01,0.02,0.01,0.01,0.001,0.001,0.948},
    mapping={{1,7},{2,2},{3,3},{4,4},{5,5},{6,6},{7,1}},
    redeclare package MediumVapour =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_SO2_HCl_HF_N2,
    redeclare record Geometry =
        ThermalSeparation.Geometry.StructuredPackedColumn.Geometry (
                                                          H=10, d=5),
    p_v_start_inlet=156000,
    p_v_start_outlet=152000,
    n_elements=3)                    annotation (Placement(transformation(extent={{-32,-14},
            {16,32}}, rotation=0)));

//    p_v_start={1.6e5,1.55e5},
  ThermalSeparation.Components.SourcesSinks.SourceGas
                                             sourceGas_Vdot(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_SO2_HCl_HF_N2,
    x={0.05,0.05,0.05,0.05,0.005,0.005,0.79},
    T=417.15)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
        origin={-24,-36})));

  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    offset=120,
    startTime=100,
    height=0) annotation (Placement(transformation(extent={{-68,-64},{-48,-44}},
          rotation=0)));
  ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink
    annotation (Placement(transformation(extent={{22,6},{36,14}}, rotation=0)));
  ThermalSeparation.Components.SourcesSinks.SplitLiquid_1p
                                                  splitLiquid_x(
 redeclare package Medium =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={24,-22})));
  ThermalSeparation.Components.SourcesSinks.CombLiquid_x
                                                combLiquid_x(
    redeclare package Medium =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O,
    x_start={0,0,0,0,0,0,1})
    annotation (Placement(transformation(
        origin={30,36},
        extent={{-10,-10},{10,10}},
        rotation=180)));

  inner SystemTS systemTS
    annotation (Placement(transformation(extent={{-64,2},{-44,22}})));
  Modelica.Blocks.Sources.Constant const(k=0.5)
    annotation (Placement(transformation(extent={{-14,-58},{6,-38}})));
equation
  connect(ramp.y, sourceGas_Vdot.Flow_in) annotation (Line(
      points={{-47,-54},{-47,-41},{-28,-41},{-28,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const.y, splitLiquid_x.split) annotation (Line(
      points={{7,-48},{12,-48},{12,-22},{18,-22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(baseStage_GF1.upStreamOut, sinkGas.gasPortIn) annotation (Line(
      points={{-24.8,29.7},{-24.8,35.85},{-24,35.85},{-24,42.4}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(baseStage_GF1.upStreamIn, sourceGas_Vdot.gasPortOut) annotation (Line(
      points={{-24.8,-11.7},{-24.8,-17.85},{-24,-17.85},{-24,-24.6}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(splitLiquid_x.liquidPortIn, baseStage_GF1.downStreamOut) annotation (
      Line(
      points={{24,-12},{16,-12},{16,-11.7},{8.8,-11.7}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(splitLiquid_x.liquidPortOut2, sinkLiquid.liquidPortIn) annotation (
      Line(
      points={{19,-32},{20,-32},{20,-52.4}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(splitLiquid_x.liquidPortOut1, combLiquid_x.liquidPortIn1) annotation (
     Line(
      points={{29,-32},{34,-32},{34,-30},{40,-30},{40,31}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(baseStage_GF1.heatPort, ambientHeatSink.heatPort1) annotation (Line(
      points={{11.68,9},{16.84,9},{16.84,10},{20.88,10}},
      color={188,51,69},
      smooth=Smooth.None));
  connect(baseStage_GF1.downStreamIn, combLiquid_x.liquidPortOut) annotation (
      Line(
      points={{8.8,29.7},{8.8,36},{20,36}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(combLiquid_x.liquidPortIn2, sourceLiquid.liquidPortOut) annotation (
      Line(
      points={{40,41},{40,56},{23.4,56}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                    graphics),
    experiment(StopTime=100, Algorithm="Dassl"),
    experimentSetupOutput,
    Documentation(info="<html>
<p><ul>
<li> spray column</li>
<li> 7 components in vapour phase </li>
<li> 7 components in liquid phase</li>
<li> absorption process</li>
</ul></p>
</html>"));
end Absorption_IdealGases_cycle_7components;
