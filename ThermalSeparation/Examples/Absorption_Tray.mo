within ThermalSeparation.Examples;
model Absorption_Tray "absorption in a tray column"
      //Lauftest: 1.2.
//Lauftest: 15.10.2011, in 15 s (laptop)

ThermalSeparation.Components.SourcesSinks.SourceLiquid
                                              sourceLiquid(
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O,
    x={0,0,1},
    T=323.15,
    Flow=0.01)             annotation (Placement(transformation(extent={{2,46},
            {22,66}}, rotation=0)));
ThermalSeparation.Components.SourcesSinks.SinkLiquid
                                            sinkLiquid(         redeclare package Medium =
                       ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O)
                      annotation (Placement(transformation(extent={{-10,-10},{
            10,10}},
                   rotation=270,
        origin={20,-62})));

ThermalSeparation.Components.SourcesSinks.SinkGas
                                         sinkGas(         redeclare package Medium =
               ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2, p=149910)
                annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-30,50})));

ThermalSeparation.Components.Columns.TrayColumn column(
    hasVapourFeed=false,
    mapping={{1,1},{2,3},{3,2}},
    redeclare package MediumVapour =
        ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2,
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O,
    nS=3,
    redeclare model Reaction = ThermalSeparation.Reaction.NoReaction,
    x_l_start_const={1e-5,1e-5,1 - 2e-5},
    x_v_start_const={0.97,0.01,0.02},
    hasLiquidFeed=false,
    entrainment=false,
    redeclare model InitOption =
        ThermalSeparation.Components.Columns.BaseClasses.Initialization.Init_x_T_p,
    h_start=0.03*ones(column.n_trays),
    redeclare record Geometry = Geometry.PlateColumn.Geometry (
        h_w=0.05,
        d=2,
        H=40),
    p_v_start_inlet=170000,
    p_v_start_outlet=150000,
    n_trays=10)          annotation (Placement(transformation(extent={{-36,-16},{12,30}},
                      rotation=0)));

  ThermalSeparation.Components.SourcesSinks.SourceGas
                                             sourceGas_Vdot(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2,
    T=363.15,
    x={0.85,0.01,0.14})
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
        origin={-28,-38})));

  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    startTime=100,
    height=0,
    offset=2.5/2)
              annotation (Placement(transformation(extent={{-84,-62},{-64,-42}},
          rotation=0)));
  ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink
    annotation (Placement(transformation(extent={{18,4},{32,14}}, rotation=0)));
  ThermalSeparation.Components.SourcesSinks.SplitLiquid_1p
                                                  splitLiquid_x(
     redeclare package Medium =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O)
                annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={24,-22})));
  ThermalSeparation.Components.SourcesSinks.CombLiquid_x
                                                combLiquid_x(
    redeclare package Medium =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_H2O,
    x_start={0,0,1})
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
      points={{-63,-52},{-52,-52},{-52,-48},{-32,-48}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const.y, splitLiquid_x.split) annotation (Line(
      points={{7,-48},{12,-48},{12,-22},{18,-22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sinkLiquid.liquidPortIn, splitLiquid_x.liquidPortOut2) annotation (
      Line(
      points={{20,-52.4},{20,-42},{20,-32},{19,-32}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(column.downStreamOut, splitLiquid_x.liquidPortIn) annotation (Line(
      points={{4.8,-13.7},{13.4,-13.7},{13.4,-12},{24,-12}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(splitLiquid_x.liquidPortOut1, combLiquid_x.liquidPortIn1) annotation (
     Line(
      points={{29,-32},{40,-32},{40,31}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(sourceLiquid.liquidPortOut, combLiquid_x.liquidPortIn2) annotation (
      Line(
      points={{23.4,56},{40,56},{40,41}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(column.downStreamIn, combLiquid_x.liquidPortOut) annotation (Line(
      points={{4.8,27.7},{14,27.7},{14,36},{20,36}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(column.heatPort, ambientHeatSink.heatPort1) annotation (Line(
      points={{8.16,7},{12,7},{12,9},{16.88,9}},
      color={188,51,69},
      smooth=Smooth.None));
  connect(column.upStreamIn, sourceGas_Vdot.gasPortOut) annotation (Line(
      points={{-28.8,-13.7},{-28.8,-19.85},{-28,-19.85},{-28,-26.6}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(column.upStreamOut, sinkGas.gasPortIn) annotation (Line(
      points={{-28.8,27.7},{-28.8,34.85},{-30,34.85},{-30,40.4}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                    graphics),
    experiment(StopTime=300, Algorithm="Dassl"),
    experimentSetupOutput,
    Documentation(info="<html>
<p><ul>
<li> tray column</li>
<li> 3 components</li>
</ul></p>
</html>"));
end Absorption_Tray;
