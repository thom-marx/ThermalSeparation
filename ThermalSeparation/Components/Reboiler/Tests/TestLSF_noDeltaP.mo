within ThermalSeparation.Components.Reboiler.Tests;
model TestLSF_noDeltaP
  // Lauftest: 25.7.2011
  //Lauftest: 14.7.2011
//Lauftest: 1.7.2011
  inner SystemTS systemTS
    annotation (Placement(transformation(extent={{60,60},{80,80}})));
  LSF_noDeltaP leanSolventFlash(
    x_v_start={0.5,0.5},
    redeclare package MediumVapour =
        ThermalSeparation.Media.IdealGasMixtures.H2O_CO2,
    inertVapour={false,false},
    nV=0,
    mapping={{1,2},{2,1}},
    length_HX=10,
    d_HX=2.9,
    redeclare model Reaction =
        ApplicationsThermalSeparation.Obsolete.CO2_Siemens (
          moleBalance=true),
    initEQ=false,
    p_start=200000)
    annotation (Placement(transformation(extent={{-38,-22},{-12,4}})));

  SourcesSinks.SinkLiquid sinkLiquid( use_Vdot=true,
    redeclare package Medium =
        ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS)
    annotation (Placement(transformation(extent={{28,16},{48,36}})));
  SourcesSinks.SinkGas sinkGas(
    use_p=true,
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.H2O_CO2,
    p=96000)
    annotation (Placement(transformation(extent={{-60,18},{-40,38}})));

  Modelica.Blocks.Sources.Ramp ramp1(
    height=0,
    offset=0,
    duration=100,
    startTime=0) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={60,-10})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow2
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={24,-10})));
  Modelica.Blocks.Sources.Ramp ramp2(
    height=0,
    duration=100,
    startTime=0,
    offset=1e-4) annotation (Placement(transformation(
        extent={{-7,-7},{7,7}},
        rotation=0,
        origin={5,19})));
ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid2(
    redeclare package MediumLiquid =
        ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS,
    use_x_in=false,
    x={0.023,0.946,0.031},
    Flow=1877*2/3600,
    T=399.15)                 annotation (Placement(transformation(extent={{-8,-50},
            {12,-30}},                                                                          rotation=0)));

equation
  connect(sinkLiquid.liquidPortIn, leanSolventFlash.liquidOut) annotation (Line(
      points={{38,33.2},{-12,33.2},{-12,40},{-23.96,40},{-23.96,1.4}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  connect(sinkGas.gasPortIn, leanSolventFlash.vapourOut) annotation (Line(
      points={{-50,21},{-50,12},{-26.82,12},{-26.82,1.4}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  connect(ramp1.y, prescribedHeatFlow2.Q_flow) annotation (Line(
      points={{49,-10},{45.25,-10},{45.25,-10},{41.5,-10},{41.5,-10},{34,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(prescribedHeatFlow2.port, leanSolventFlash.heatPort) annotation (Line(
      points={{14,-10},{0,-10},{0,-10.82},{-14.6,-10.82}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(ramp2.y, sinkLiquid.Vdot_set) annotation (Line(
      points={{12.7,19},{24,19},{24,29.3},{30.9,29.3}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sourceLiquid2.liquidPortOut, leanSolventFlash.liquidIn) annotation (
     Line(
      points={{2,-48.6},{2,-60},{-25,-60},{-25,-19.4}},
      color={0,0,0},
      thickness=1,
      smooth=Smooth.None));
  annotation (Diagram(graphics));
end TestLSF_noDeltaP;
