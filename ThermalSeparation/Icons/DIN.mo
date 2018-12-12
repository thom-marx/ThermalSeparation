within ThermalSeparation.Icons;
package DIN
  model HeatExchanger

    LiquidPortOut liquidPortOut
      annotation (Placement(transformation(extent={{-100,-10},{-80,10}}),
          iconTransformation(extent={{-100,-10},{-80,10}})));
    LiquidPortIn liquidPortIn
      annotation (Placement(transformation(extent={{80,-10},{100,10}}),
          iconTransformation(extent={{80,-10},{100,10}})));
    LiquidPortIn liquidPortIn1
      annotation (Placement(transformation(extent={{-60,70},{-40,90}}),
          iconTransformation(extent={{-60,70},{-40,90}})));
    LiquidPortOut liquidPortOut1
      annotation (Placement(transformation(extent={{40,70},{60,90}}),
          iconTransformation(extent={{40,70},{60,90}})));
  equation
    connect(liquidPortOut, liquidPortOut) annotation (Line(
        points={{-90,0},{-56,0},{-56,0},{-90,0}},
        color={0,0,0},
        smooth=Smooth.None));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{
              -100,-100},{100,100}}),
                        graphics),                                  Icon(
          coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}),
          graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={0,0,0}),
            Line(
            points={{-40,80},{-40,-40},{0,40},{40,-40},{40,80}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end HeatExchanger;

  model Compressors

    GasPortIn gasPortIn annotation (Placement(transformation(extent={{-10,
              -100},{10,-80}}),
                          iconTransformation(extent={{-10,-100},{10,-80}})));
    GasPortOut gasPortOut annotation (Placement(transformation(extent={{-10,80},
              {10,100}}),    iconTransformation(extent={{-10,80},{10,100}})));
  equation
    connect(gasPortOut, gasPortOut) annotation (Line(
        points={{0,90},{0,90},{0,90},{0,90}},
        color={255,255,0},
        smooth=Smooth.None));
    annotation (Diagram(graphics), Icon(graphics={
          Ellipse(
            extent={{-80,-80},{80,80}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-32,74},{-76,-28}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{32,74},{76,-28}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end Compressors;

  model Pumps

    LiquidPortIn liquidPortIn annotation (Placement(transformation(extent={{-10,80},
              {10,100}}),        iconTransformation(extent={{-10,80},{10,100}})));
    LiquidPortOut liquidPortOut annotation (Placement(transformation(extent={{-10,
              -100},{10,-80}}),     iconTransformation(extent={{-10,-100},{10,
              -80}})));
    annotation (Diagram(graphics), Icon(graphics={Ellipse(
            extent={{-80,-80},{80,80}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Line(
            points={{-80,0},{0,-80},{80,0}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end Pumps;

  model Tank

    LiquidPortIn liquidPortIn
      annotation (Placement(transformation(extent={{-60,110},{-40,130}}),
          iconTransformation(extent={{-60,110},{-40,130}})));
    LiquidPortIn liquidPortIn1
      annotation (Placement(transformation(extent={{120,110},{140,130}}),
          iconTransformation(extent={{120,110},{140,130}})));
    LiquidPortOut liquidPortOut
      annotation (Placement(transformation(extent={{30,-100},{50,-80}}),
          iconTransformation(extent={{30,-100},{50,-80}})));
    annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
              -100},{180,180}}),
                     graphics={Rectangle(extent={{-40,160},{120,-80}},
              lineColor={0,0,0})}), Diagram(coordinateSystem(
            preserveAspectRatio=true, extent={{-100,-100},{180,180}})));
  end Tank;

  model Flash

    LiquidPortIn liquidPortIn
      annotation (Placement(transformation(extent={{40,-80},{60,-60}}),
          iconTransformation(extent={{40,-80},{60,-60}})));
    GasPortOut gasPortOut
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
    LiquidPortOut liquidPortOut
      annotation (Placement(transformation(extent={{60,60},{80,80}}),
          iconTransformation(extent={{60,60},{80,80}})));
  equation
    connect(gasPortOut, gasPortOut) annotation (Line(
        points={{-70,70},{-70,70}},
        color={0,0,255},
        smooth=Smooth.None));
    annotation (Diagram(graphics), Icon(graphics={Polygon(
            points={{-80,80},{80,80},{40,-80},{-40,-80},{-80,80}},
            lineColor={0,0,0},
            smooth=Smooth.None)}));
  end Flash;

  model Reboiler

    LiquidPortIn liquidPortIn
      annotation (Placement(transformation(extent={{60,60},{80,80}})));
    GasPortOut gasPortOut
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
    LiquidPortOut liquidPortOut
      annotation (Placement(transformation(extent={{60,-80},{80,-60}})));
    annotation (Diagram(graphics), Icon(graphics={Polygon(
            points={{-80,-80},{-60,80},{60,80},{80,-80},{-80,-80}},
            lineColor={0,0,0},
            smooth=Smooth.None), Line(
            points={{-60,-60},{-50,-40},{-40,-60},{-30,-40},{-20,-60},{-10,
                -40},{0,-60},{10,-40},{20,-60},{30,-40},{40,-60},{50,-40},{60,
                -60}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end Reboiler;

  model Sump

    LiquidPortIn liquidPortIn2
      annotation (Placement(transformation(extent={{-60,110},{-40,130}}),
          iconTransformation(extent={{-60,110},{-40,130}})));
    LiquidPortIn liquidPortIn3
      annotation (Placement(transformation(extent={{120,110},{140,130}}),
          iconTransformation(extent={{120,110},{140,130}})));
    LiquidPortOut liquidPortOut1
      annotation (Placement(transformation(extent={{30,-100},{50,-80}}),
          iconTransformation(extent={{30,-100},{50,-80}})));
    annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
              -100},{180,180}}),
                     graphics={Rectangle(extent={{-40,160},{120,-80}},
              lineColor={0,0,0})}), Diagram(coordinateSystem(
            preserveAspectRatio=true, extent={{-100,-100},{180,180}})));
  end Sump;

  model LiquidSplitter

    LiquidPortIn liquidPortIn1
      annotation (Placement(transformation(extent={{-80,-10},{-60,10}}),
          iconTransformation(extent={{-80,-10},{-60,10}})));
    LiquidPortOut liquidPortOut
      annotation (Placement(transformation(extent={{60,-80},{80,-60}}),
          iconTransformation(extent={{60,-80},{80,-60}})));
    LiquidPortOut liquidPortOut1 annotation (Placement(transformation(extent={{60,60},
              {80,80}}),          iconTransformation(extent={{60,60},{80,80}})));
    annotation (Diagram(graphics), Icon(graphics={
          Line(
            points={{-27,63},{29,1},{-29,-63}},
            color={0,0,0},
            smooth=Smooth.None,
            origin={49,1},
            rotation=180),
          Line(
            points={{-66,33},{-4,-33},{66,-33}},
            color={0,0,0},
            smooth=Smooth.None,
            origin={-4,-43},
            rotation=180),
          Line(
            points={{-66,-33},{-4,33},{66,33}},
            color={0,0,0},
            smooth=Smooth.None,
            origin={-4,43},
            rotation=180)}));
  end LiquidSplitter;

  connector GasPortIn

    annotation (Icon(graphics={Ellipse(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Forward)}));
  end GasPortIn;

  connector GasPortOut

    annotation (Icon(graphics={Ellipse(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Forward)}));
  end GasPortOut;

  connector LiquidPortIn

    annotation (Icon(graphics={Ellipse(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Backward)}));
  end LiquidPortIn;

  connector LiquidPortOut

    annotation (Icon(graphics={Ellipse(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Backward)}));
  end LiquidPortOut;

  connector HeatPort

    annotation (Icon(graphics={Ellipse(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.CrossDiag)}));
  end HeatPort;

  model LiquidSource

    LiquidPortOut liquidPortOut
      annotation (Placement(transformation(extent={{60,60},{80,80}}),
          iconTransformation(extent={{-58,-58},{-38,-38}})));
  equation
    connect(liquidPortOut, liquidPortOut) annotation (Line(
        points={{70,70},{44,70},{44,70},{70,70}},
        color={0,0,255},
        smooth=Smooth.None));
    annotation (Diagram(graphics), Icon(graphics={
                               Ellipse(
            extent={{-60,-60},{60,60}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Backward)}));
  end LiquidSource;

  model LiquidSink

    LiquidPortIn liquidPortIn
      annotation (Placement(transformation(extent={{-198,-102},{-38,58}}),
          iconTransformation(extent={{-58,38},{-38,58}})));
    annotation (Diagram(graphics), Icon(graphics={
                               Ellipse(
            extent={{-60,-60},{60,60}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Backward)}));
  end LiquidSink;

  model GasSource

    GasPortOut gasPortOut annotation (Placement(transformation(extent={{38,38},
              {58,58}}), iconTransformation(extent={{38,38},{58,58}})));
    annotation (Diagram(graphics), Icon(graphics={
                               Ellipse(
            extent={{-60,-60},{60,60}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Forward)}));
  end GasSource;

  model GasSink

    GasPortIn gasPortIn
      annotation (Placement(transformation(extent={{38,-58},{58,-38}}),
          iconTransformation(extent={{38,-58},{58,-38}})));
    annotation (Icon(graphics={Ellipse(
            extent={{-60,-60},{60,60}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Forward)}), Diagram(graphics));
  end GasSink;

  model AmbientHeatSink

    HeatPort heatPort
      annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
    annotation (Diagram(graphics), Icon(graphics={Line(
            points={{-60,12},{-36,-64},{-20,78},{0,-80},{20,80},{40,-80},{60,
                78},{78,-64},{94,10}},
            color={0,0,0},
            smooth=Smooth.Bezier,
            arrow={Arrow.None,Arrow.Filled})}));
  end AmbientHeatSink;

  model GasSplitter

    GasPortOut gasPortOut
      annotation (Placement(transformation(extent={{-80,-80},{-60,-60}}),
          iconTransformation(extent={{-80,-80},{-60,-60}})));
    GasPortIn gasPortIn
      annotation (Placement(transformation(extent={{60,-10},{80,10}}),
          iconTransformation(extent={{60,-10},{80,10}})));
    GasPortOut gasPortOut1
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
    annotation (Diagram(graphics), Icon(graphics={
          Line(
            points={{-76,62},{-20,0},{-78,-64}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{-62,-76},{0,-10},{70,-10}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{-62,76},{0,10},{70,10}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end GasSplitter;

  model HeatExchangeH2O

    LiquidPortOut liquidPortOut
      annotation (Placement(transformation(extent={{-100,-10},{-80,10}}),
          iconTransformation(extent={{-100,-10},{-80,10}})));
    LiquidPortIn liquidPortIn
      annotation (Placement(transformation(extent={{80,-10},{100,10}}),
          iconTransformation(extent={{80,-10},{100,10}})));
    Modelica.Fluid.Interfaces.FluidPorts_b ports_b
      annotation (Placement(transformation(extent={{-80,20},{-60,100}}),
          iconTransformation(extent={{-10,-40},{10,40}},
          rotation=90,
          origin={-60,90})));
    Modelica.Fluid.Interfaces.FluidPorts_b ports_b1
      annotation (Placement(transformation(extent={{60,20},{80,100}}),
          iconTransformation(extent={{-10,-40},{10,40}},
          rotation=90,
          origin={60,90})));
  equation
    connect(ports_b, ports_b) annotation (Line(
        points={{-70,60},{-70,76},{-70,60},{-70,60}},
        color={0,127,255},
        smooth=Smooth.None));
    annotation (Icon(graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor=
                {0,0,0}), Line(
            points={{-40,80},{-40,-40},{0,40},{40,-40},{40,80}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end HeatExchangeH2O;

  model Columns

    GasPortIn gasPortIn
      annotation (Placement(transformation(extent={{30,-100},{50,-80}}),
          iconTransformation(extent={{30,-100},{50,-80}})));
    GasPortOut gasPortOut annotation (Placement(transformation(extent={{30,160},
              {50,180}}),     iconTransformation(extent={{30,160},{50,180}})));
    LiquidPortIn liquidPortIn
      annotation (Placement(transformation(extent={{120,130},{140,150}}),
          iconTransformation(extent={{120,130},{140,150}})));
    LiquidPortOut liquidPortOut
      annotation (Placement(transformation(extent={{120,-70},{140,-50}}),
          iconTransformation(extent={{120,-70},{140,-50}})));
    ThermalSeparation.Icons.DIN.HeatPort ambientHeatSink
      annotation (Placement(transformation(extent={{120,30},{140,50}}),
          iconTransformation(extent={{120,30},{140,50}})));
    annotation (Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{180,180}},
          grid={2,2},
          initialScale=0.2), graphics), Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{180,180}},
          grid={2,2},
          initialScale=0.2), graphics={
          Rectangle(extent={{-40,160},{120,-80}}, lineColor={0,0,0}),
          Line(
            points={{-40,-40},{120,-40}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{-40,120},{120,-40}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{-40,120},{120,120}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{120,120},{-40,-40}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end Columns;

  model LiquidMixer

    LiquidPortIn liquidPortIn
      annotation (Placement(transformation(extent={{-70,80},{-50,100}}),
          iconTransformation(extent={{-70,80},{-50,100}})));
    LiquidPortIn liquidPortIn1
      annotation (Placement(transformation(extent={{-30,80},{-10,100}}),
          iconTransformation(extent={{-30,80},{-10,100}})));
    LiquidPortOut liquidPortOut
      annotation (Placement(transformation(extent={{90,-60},{110,-40}}),
          iconTransformation(extent={{90,-60},{110,-40}})));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{
              -100,-100},{140,140}}),
                        graphics), Icon(coordinateSystem(preserveAspectRatio=
              true, extent={{-100,-100},{140,140}}),
                                        graphics={
          Rectangle(extent={{-80,80},{120,-40}}, lineColor={0,0,0}),
          Line(
            points={{-100,20},{120,20}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{-60,40},{-40,60},{0,-20},{20,0}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{20,40},{40,60},{80,-20},{100,0}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end LiquidMixer;

  model GasMixer

    GasPortOut gasPortOut
      annotation (Placement(transformation(extent={{-80,-10},{-60,10}}),
          iconTransformation(extent={{-80,-10},{-60,10}})));
    GasPortIn gasPortIn
      annotation (Placement(transformation(extent={{60,-80},{80,-60}}),
          iconTransformation(extent={{60,-80},{80,-60}})));
    GasPortIn gasPortIn1
      annotation (Placement(transformation(extent={{60,60},{80,80}}),
          iconTransformation(extent={{60,60},{80,80}})));
    annotation (Diagram(graphics), Icon(graphics={
          Line(
            points={{-27,63},{29,1},{-29,-63}},
            color={0,0,0},
            smooth=Smooth.None,
            origin={49,1},
            rotation=180),
          Line(
            points={{-66,-33},{-4,33},{66,33}},
            color={0,0,0},
            smooth=Smooth.None,
            origin={-4,43},
            rotation=180),
          Line(
            points={{-66,33},{-4,-33},{66,-33}},
            color={0,0,0},
            smooth=Smooth.None,
            origin={-4,-43},
            rotation=180)}));
  end GasMixer;

  model FGCooler2

    ThermalSeparation.Icons.DIN.Pumps pumps
      annotation (Placement(transformation(extent={{-20,54},{0,74}})));
    ThermalSeparation.Icons.DIN.Compressors compressors
      annotation (Placement(transformation(extent={{-50,20},{-30,40}})));
    ThermalSeparation.Icons.DIN.Tank tank
      annotation (Placement(transformation(extent={{-20,80},{0,100}})));
    ThermalSeparation.Icons.DIN.LiquidSink liquidSink
      annotation (Placement(transformation(extent={{-20,-80},{0,-60}})));
    ThermalSeparation.Icons.DIN.GasSource gasSource
      annotation (Placement(transformation(extent={{-80,-80},{-60,-60}})));
    ThermalSeparation.Icons.DIN.GasSink gasSink
      annotation (Placement(transformation(extent={{-62,140},{-42,160}})));
    ThermalSeparation.Icons.DIN.LiquidSource liquidSource
      annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
    ThermalSeparation.Icons.DIN.LiquidSource liquidSource1
      annotation (Placement(transformation(extent={{0,90},{20,110}})));
    ThermalSeparation.Icons.DIN.LiquidSink liquidSink1
      annotation (Placement(transformation(extent={{80,-30},{100,-10}})));
    ThermalSeparation.Icons.DIN.LiquidSplitter splitter1
      annotation (Placement(transformation(extent={{100,38},{120,58}})));
    ThermalSeparation.Icons.DIN.HeatExchanger heatExchanger1
      annotation (Placement(transformation(extent={{140,-10},{160,10}})));
    ThermalSeparation.Icons.DIN.Pumps pumps1
      annotation (Placement(transformation(extent={{224,0},{244,-20}})));
    ThermalSeparation.Icons.DIN.Compressors compressors1
      annotation (Placement(transformation(extent={{294,-30},{314,-10}})));
    ThermalSeparation.Icons.DIN.Flash flash
      annotation (Placement(transformation(extent={{300,-60},{280,-40}})));
    ThermalSeparation.Icons.DIN.Pumps pumps2
      annotation (Placement(transformation(extent={{320,-60},{340,-40}})));
    ThermalSeparation.Icons.DIN.Sump sump
      annotation (Placement(transformation(extent={{320,0},{340,20}})));
    ThermalSeparation.Icons.DIN.Reboiler reboiler
      annotation (Placement(transformation(extent={{320,-20},{340,0}})));
    ThermalSeparation.Icons.DIN.LiquidSource liquidSource2
      annotation (Placement(transformation(extent={{320,164},{340,184}})));
    ThermalSeparation.Icons.DIN.GasSink gasSink1
      annotation (Placement(transformation(extent={{280,140},{300,160}})));
    ThermalSeparation.Icons.DIN.AmbientHeatSink ambientHeatSink
      annotation (Placement(transformation(extent={{-20,-48},{0,-28}})));
    ThermalSeparation.Icons.DIN.AmbientHeatSink ambientHeatSink2
      annotation (Placement(transformation(extent={{320,110},{340,130}})));
    ThermalSeparation.Icons.DIN.AmbientHeatSink ambientHeatSink3
      annotation (Placement(transformation(extent={{320,30},{340,50}})));
    ThermalSeparation.Icons.DIN.Columns columns
      annotation (Placement(transformation(extent={{-60,100},{-20,140}})));
    ThermalSeparation.Icons.DIN.Columns columns1
      annotation (Placement(transformation(extent={{-60,-60},{-20,-20}})));
    ThermalSeparation.Icons.DIN.Columns columns2
      annotation (Placement(transformation(extent={{280,100},{320,140}})));
    ThermalSeparation.Icons.DIN.Columns columns3
      annotation (Placement(transformation(extent={{280,20},{320,60}})));
    Modelica.Fluid.Sources.MassFlowSource_h boundary(nPorts=1)
      annotation (Placement(transformation(extent={{0,160},{20,180}})));
    Modelica.Fluid.Sources.FixedBoundary boundary1(nPorts=1)
      annotation (Placement(transformation(extent={{60,180},{80,200}})));
    ThermalSeparation.Icons.DIN.HeatExchangeH2O heatExchangeH2O
      annotation (Placement(transformation(extent={{40,140},{60,160}})));
    ThermalSeparation.Icons.DIN.AmbientHeatSink ambientHeatSink1
      annotation (Placement(transformation(extent={{-20,112},{0,132}})));
    ThermalSeparation.Icons.DIN.LiquidSplitter liquidSplitter
      annotation (Placement(transformation(extent={{102,-10},{82,10}})));
    ThermalSeparation.Icons.DIN.LiquidMixer liquidMixer
      annotation (Placement(transformation(extent={{308,78},{328,98}})));
    ThermalSeparation.Icons.DIN.LiquidMixer liquidMixer1
      annotation (Placement(transformation(extent={{320,142},{340,162}})));
    ThermalSeparation.Icons.DIN.GasMixer gasMixer
      annotation (Placement(transformation(extent={{286,-6},{306,14}})));
  equation
    connect(splitter1.liquidPortOut, heatExchanger1.liquidPortIn1)
      annotation (Line(
        points={{117,41},{145.5,41},{145.5,8},{145,8}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(compressors1.gasPortIn, flash.gasPortOut) annotation (Line(
        points={{304,-29},{306.5,-29},{306.5,-43},{297,-43}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(pumps.liquidPortOut, splitter1.liquidPortIn1) annotation (Line(
        points={{-10,55},{-10,48.5},{103,48.5},{103,48}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(columns3.liquidPortOut, sump.liquidPortIn1) annotation (Line(
        points={{314,26},{320,26},{320,17},{323,17}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(columns1.gasPortIn, gasSource.gasPortOut) annotation (Line(
        points={{-40,-58.5714},{-66,-58.5714},{-66,-65.2},{-65.2,-65.2}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(compressors.gasPortOut, columns.gasPortIn) annotation (Line(
        points={{-40,39},{-40,101.429}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(compressors.gasPortIn, columns1.gasPortOut) annotation (Line(
        points={{-40,21},{-40,-21.4286}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(heatExchangeH2O.ports_b, boundary.ports[1]) annotation (Line(
        points={{44,159},{44,162.5},{20,162.5},{20,170}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(boundary1.ports[1], heatExchangeH2O.ports_b1) annotation (Line(
        points={{80,190},{70,190},{70,159},{56,159}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(heatExchangeH2O.liquidPortOut, columns.liquidPortIn) annotation (
        Line(
        points={{41,150},{-11.5,150},{-11.5,134.286},{-27.1429,134.286}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(ambientHeatSink1.heatPort, columns.ambientHeatSink) annotation (
        Line(
        points={{-17,123},{-21.5,123},{-21.5,120},{-27.1429,120}},
        color={0,0,0},
        smooth=Smooth.Bezier));
    connect(columns1.ambientHeatSink, ambientHeatSink.heatPort) annotation (
        Line(
        points={{-27.1429,-40},{-17,-40},{-17,-37}},
        color={0,0,0},
        smooth=Smooth.Bezier));
    connect(ambientHeatSink3.heatPort, columns3.ambientHeatSink) annotation (
        Line(
        points={{323,41},{318.5,41},{318.5,40},{312.857,40}},
        color={0,0,0},
        smooth=Smooth.Bezier));
    connect(ambientHeatSink2.heatPort, columns2.ambientHeatSink) annotation (
        Line(
        points={{323,121},{319.5,121},{319.5,120},{312.857,120}},
        color={0,0,0},
        smooth=Smooth.Bezier));
    connect(gasMixer.gasPortIn, compressors1.gasPortOut) annotation (Line(
        points={{303,-3},{303,-9.5},{304,-9.5},{304,-11}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(liquidSplitter.liquidPortIn1, heatExchanger1.liquidPortOut)
      annotation (Line(
        points={{99,0},{141,0}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(heatExchangeH2O.liquidPortIn, liquidSplitter.liquidPortOut1)
      annotation (Line(
        points={{59,150},{68,150},{68,7},{85,7}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(tank.liquidPortIn1, liquidSource1.liquidPortOut) annotation (Line(
        points={{-3.57143,95.7143},{0.5,95.7143},{0.5,95.2},{5.2,95.2}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(liquidSink.liquidPortIn, columns1.liquidPortOut) annotation (Line(
        points={{-14.8,-65.2},{-14.8,-59.6},{-27.1429,-59.6},{-27.1429,
            -54.2857}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(gasMixer.gasPortIn1, reboiler.gasPortOut) annotation (Line(
        points={{303,11},{310.5,11},{310.5,-3},{323,-3}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(reboiler.liquidPortIn, sump.liquidPortOut1) annotation (Line(
        points={{337,-3},{337,-0.5},{337,0.714286},{330,0.714286}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(pumps2.liquidPortIn, reboiler.liquidPortOut) annotation (Line(
        points={{330,-41},{338,-41},{338,-18},{337,-17}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(columns3.liquidPortIn, liquidMixer.liquidPortOut) annotation (
        Line(
        points={{312.857,54.2857},{346,54.2857},{346,82.1667},{324.667,
            82.1667}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(heatExchanger1.liquidPortIn, pumps1.liquidPortOut) annotation (
        Line(
        points={{159,0},{223.5,0},{223.5,-1},{234,-1}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(flash.liquidPortIn, pumps2.liquidPortOut) annotation (Line(
        points={{285,-57},{312.5,-57},{312.5,-59},{330,-59}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(pumps1.liquidPortIn, flash.liquidPortOut) annotation (Line(
        points={{234,-19},{254,-19},{254,-43},{283,-43}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(sump.liquidPortIn3, reboiler.liquidPortOut) annotation (Line(
        points={{336.429,15.7143},{336.429,16},{368,16},{368,-17},{337,-17}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(heatExchanger1.liquidPortOut1, liquidMixer.liquidPortIn)
      annotation (Line(
        points={{155,8},{272,8},{272,93.8333},{311.333,93.8333}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(columns2.liquidPortOut, liquidMixer.liquidPortIn1) annotation (
        Line(
        points={{312.857,105.714},{314.667,105.714},{314.667,93.8333}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(columns2.gasPortIn, columns3.gasPortOut) annotation (Line(
        points={{300,101.429},{300,58.5714}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(columns1.liquidPortIn, liquidSource.liquidPortOut) annotation (
        Line(
        points={{-27.1429,-25.7143},{-21.5714,-25.7143},{-21.5714,-14.8},{
            -14.8,-14.8}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(columns.liquidPortOut, tank.liquidPortIn) annotation (Line(
        points={{-27.1429,105.714},{-21.5714,105.714},{-21.5714,95.7143},{
            -16.4286,95.7143}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(tank.liquidPortOut, pumps.liquidPortIn) annotation (Line(
        points={{-10,80.7143},{-10,73}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(liquidSplitter.liquidPortOut, liquidSink1.liquidPortIn)
      annotation (Line(
        points={{85,-7},{85,-11.5},{85.2,-11.5},{85.2,-15.2}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(splitter1.liquidPortOut1, liquidMixer1.liquidPortIn) annotation (
        Line(
        points={{117,55},{144.5,55},{144.5,157.833},{323.333,157.833}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(gasSink1.gasPortIn, columns2.gasPortOut) annotation (Line(
        points={{294.8,145.2},{298.4,141.6},{298.4,138.571},{300,138.571}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(liquidMixer1.liquidPortIn1, liquidSource2.liquidPortOut)
      annotation (Line(
        points={{326.667,157.833},{326.667,162.916},{325.2,162.916},{325.2,
            169.2}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(liquidMixer1.liquidPortOut, columns2.liquidPortIn) annotation (
        Line(
        points={{336.667,146.167},{325.333,146.167},{325.333,134.286},{
            312.857,134.286}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(gasSink.gasPortIn, columns.gasPortOut) annotation (Line(
        points={{-47.2,145.2},{-44.6,145.2},{-44.6,138.571},{-40,138.571}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(gasMixer.gasPortOut, columns3.gasPortIn) annotation (Line(
        points={{289,4},{288.5,4},{288.5,21.4286},{300,21.4286}},
        color={0,0,0},
        smooth=Smooth.None));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
              -100},{400,200}}),      graphics), Icon(coordinateSystem(
            preserveAspectRatio=true, extent={{-100,-100},{400,200}})));
  end FGCooler2;
end DIN;
