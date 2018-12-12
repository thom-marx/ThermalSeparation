within ThermalSeparation.Components.Reboiler.Tests;
model TestLSF
// Lauftest: 25.7.2011
// Lauftest: 14.7.2011
  inner SystemTS systemTS
    annotation (Placement(transformation(extent={{60,60},{80,80}})));
  LSF_deltaP leanSolventFlash(
    x_v_start={0.5,0.5},
    redeclare package MediumVapour =
        ThermalSeparation.Media.IdealGasMixtures.H2O_CO2,
    inertVapour={false,false},
    nV=0,
    mapping={{1,2},{2,1}},
    length_HX=10,
    d_HX=2.9,
    initEQ=true,
    redeclare package MediumLiquid =
        ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS,
    p_start=200000,
    redeclare model Reaction =
        ApplicationsThermalSeparation.Reaction.CO2_Siemens)
    annotation (Placement(transformation(extent={{-54,-22},{-28,4}})));

  SourcesSinks.SinkLiquid sinkLiquid(use_Vdot=true,
    redeclare package Medium =
        ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS,
    use_p=false)
    annotation (Placement(transformation(extent={{12,16},{32,36}})));
  SourcesSinks.SinkGas sinkGas(
    use_p=true,
    p=96000)
    annotation (Placement(transformation(extent={{-76,18},{-56,38}})));
  Modelica.Blocks.Sources.Ramp ramp1(
    height=0,
    offset=0,
    duration=100,
    startTime=0) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={44,-10})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow2
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={8,-10})));
  Modelica.Blocks.Sources.Ramp ramp2(
    height=0,
    duration=100,
    startTime=0,
    offset=1e-4) annotation (Placement(transformation(
        extent={{-7,-7},{7,7}},
        rotation=0,
        origin={-11,19})));
ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid2(
    redeclare package MediumLiquid =
        ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS,
    use_x_in=false,
    x={0.023,0.946,0.031},
    T=399.15,
    Flow=1877/3600*2)         annotation (Placement(transformation(extent={{-24,-50},
            {-4,-30}},                                                                          rotation=0)));
equation

  connect(sinkLiquid.liquidPortIn, leanSolventFlash.liquidOut) annotation (Line(
      points={{22,33.2},{-28,33.2},{-28,40},{-39.96,40},{-39.96,1.4}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  connect(sinkGas.gasPortIn, leanSolventFlash.vapourOut) annotation (Line(
      points={{-66,21},{-66,12},{-42.82,12},{-42.82,1.4}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  connect(ramp1.y, prescribedHeatFlow2.Q_flow) annotation (Line(
      points={{33,-10},{29.25,-10},{29.25,-10},{25.5,-10},{25.5,-10},{18,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(prescribedHeatFlow2.port, leanSolventFlash.heatPort) annotation (Line(
      points={{-2,-10},{-16,-10},{-16,-10.82},{-30.6,-10.82}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(ramp2.y, sinkLiquid.Vdot_set) annotation (Line(
      points={{-3.3,19},{8,19},{8,29.3},{14.9,29.3}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sourceLiquid2.liquidPortOut, leanSolventFlash.liquidIn) annotation (
     Line(
      points={{-14,-48.6},{-14,-60},{-41,-60},{-41,-19.4}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  annotation (Diagram(graphics));
end TestLSF;
