within ThermalSeparation.Components.Compressors.Tests;
model CompressorSimple
  import ThermalSeparation;
  // Lauftest: 14.7.2011
    //X_start_in=reference_X)
  ThermalSeparation.Components.Compressors.CompressorSimple
                          Compressor(redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.CO2_H2O,
    P_drive_const=10000,
    useP=true)
    annotation (Placement(transformation(extent={{-30,-6},{-10,14}})));
   // use_alpha_in=false,

  SourcesSinks.SinkGas sinkGas(redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.CO2_H2O, p=205000)
    annotation (Placement(transformation(extent={{18,20},{38,40}})));
  SourcesSinks.SourceGas sourceGas_Vdot(
    x={0.9,0.1},
    redeclare package Medium = ThermalSeparation.Media.IdealGasMixtures.CO2_H2O,
    use_Flow=false,
    flowOption=ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowVdot,
    T=323.15,
    Flow=0.55)
              annotation (Placement(transformation(extent={{-74,6},{-54,26}})));

  Modelica.Blocks.Sources.Ramp ramp(
    height=40000,
    duration=0.5,
    offset=10000,
    startTime=0.2)
    annotation (Placement(transformation(extent={{-76,-24},{-56,-4}})));
equation
  connect(Compressor.gasPortIn, sourceGas_Vdot.gasPortOut) annotation (Line(
      points={{-20,-5.4},{-20,34},{-52.6,34},{-52.6,16}},
      color={160,160,0},
      thickness=1,
      smooth=Smooth.None));
  connect(Compressor.gasPortOut, sinkGas.gasPortIn) annotation (Line(
      points={{-20,13.4},{-20,-22},{18.4,-22},{18.4,30}},
      color={160,160,0},
      thickness=1,
      smooth=Smooth.None));
  connect(ramp.y, Compressor.u) annotation (Line(
      points={{-55,-14},{-42,-14},{-42,4},{-32,4}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},
            {100,100}}), graphics), Commands);
end CompressorSimple;
