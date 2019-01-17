within ThermalSeparation.Components.Columns.BaseClasses;
model InternalStepSequence_StateGraph1
  Modelica.StateGraph.InitialStep ColdState(
    nOut=1,
    nIn=1)
    annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=270,
        origin={-32,36})));
  Modelica.StateGraph.TransitionWithSignal T1
    annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=270,
        origin={-32,12})));
  Modelica.StateGraph.Step HeatUpState(
    nOut=1,
    nIn=1) annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=270,
        origin={-32,-14})));
 Modelica.StateGraph.TransitionWithSignal T2
    annotation (Placement(transformation(extent={{-4,-4},{4,4}},
        rotation=270,
        origin={-32,-38})));
  Modelica.StateGraph.Step NormalOperation(
    nOut=1,
    nIn=1) annotation (Placement(transformation(
        extent={{-4,-4},{4,4}},
        rotation=0,
        origin={-4,-54})));
  Modelica.StateGraph.TransitionWithSignal T3(
    waitTime=1)                                              annotation (Placement(
        transformation(
        extent={{-4,-4},{4,4}},
        rotation=90,
        origin={20,-40})));
  Modelica.StateGraph.Step ShutDown(
    nIn=1,
    nOut=1) annotation (Placement(transformation(
        extent={{-4,-4},{4,4}},
        rotation=90,
        origin={20,12})));
  Modelica.StateGraph.TransitionWithSignal T4 annotation (Placement(
        transformation(
        extent={{-4,-4},{4,4}},
        rotation=180,
        origin={-2,60})));
  Modelica.Blocks.Interfaces.BooleanInput shutDown_signal annotation (Placement(
        transformation(
        extent={{-16,-16},{16,16}},
        rotation=90,
        origin={38,-100})));
  Modelica.Blocks.Interfaces.BooleanInput StartUp_signal1 annotation (Placement(
        transformation(
        extent={{-17,-17},{17,17}},
        rotation=90,
        origin={81,-99})));
  Modelica.Blocks.Interfaces.BooleanInput endPhase1 annotation (Placement(
      transformation(
      extent={{-16,-16},{16,16}},
      rotation=0,
      origin={-98,12})));
  Modelica.Blocks.Interfaces.BooleanInput endPhase2 annotation (Placement(
      transformation(
      extent={{-2,-2},{2,2}},
      rotation=0,
      origin={-44,-38})));
      inner Modelica.StateGraph.StateGraphRoot stateGraphRoot;
equation
  connect(ColdState.outPort[1], T1.inPort)
    annotation (Line(points={{-32,31.8},{-32,26},{-32,13.6},{-32,13.6}},
                                                   color={0,0,0}));
  connect(HeatUpState.inPort[1], T1.outPort)
    annotation (Line(points={{-32,-9.6},{-32,-9.6},{-32,11.4}},
                                                          color={0,0,0}));
  connect(T2.inPort, HeatUpState.outPort[1])
    annotation (Line(points={{-32,-36.4},{-32,-18.2},{-32,-18.2}},
                                                     color={0,0,0}));
  connect(NormalOperation.inPort[1], T2.outPort)
    annotation (Line(points={{-8.4,-54},{-32,-54},{-32,-38.6}},
                                                            color={0,0,0}));
  connect(NormalOperation.outPort[1], T3.inPort) annotation (Line(points={{0.2,-54},{0.2,-54},{20,-54},{20,-41.6}},
                                        color={0,0,0}));
  connect(ShutDown.inPort[1], T3.outPort)
    annotation (Line(points={{20,7.6},{20,7.6},{20,-39.4}},
                                                        color={0,0,0}));
  connect(T4.outPort, ColdState.inPort[1])
    annotation (Line(points={{-2.6,60},{-32,60},{-32,40.4}},
                                                         color={0,0,0}));
  connect(T4.inPort, ShutDown.outPort[1])
    annotation (Line(points={{-0.4,60},{20,60},{20,16.2}},
                                                        color={0,0,0}));
  connect(endPhase1, T1.condition) annotation (Line(points={{-98,12},{-36.8,12}}, color={255,0,255}));
  connect(endPhase2, T2.condition) annotation (Line(points={{-44,-38},{-36.8,-38}}, color={255,0,255}));
  connect(shutDown_signal, T3.condition) annotation (Line(points={{38,-100},{38,-40},{24.8,-40}}, color={255,0,255}));
  connect(StartUp_signal1, T4.condition) annotation (Line(points={{81,-99},{81,82},{-2,82},{-2,64.8}}, color={255,0,255}));
                                           annotation (Placement(transformation(extent={{-80,60},{-60,80}})),
              Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-94,94},{94,-94}},
          lineColor={85,170,255},
          lineThickness=1),
        Rectangle(
          extent={{-12,60},{12,40}},
          lineColor={135,135,135},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          radius=60,
          lineThickness=0.5),
        Rectangle(
          extent={{-12,0},{12,-20}},
          lineColor={135,135,135},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          radius=60,
          lineThickness=0.5),
        Rectangle(
          extent={{-12,-60},{12,-80}},
          lineColor={135,135,135},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          radius=60,
          lineThickness=0.5),
        Rectangle(
          extent={{-30,19},{32,24}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          radius=10),
        Rectangle(
          extent={{-30,-41},{32,-36}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          radius=10),
        Line(
          points={{0,40},{0,24}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{0,20},{0,0}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{0,-20},{0,-36},{0,-38}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{0,-40},{0,-56},{0,-60}},
          color={0,0,0},
          thickness=0.5),
        Ellipse(
          extent={{-6,-58},{6,-62}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-6,2},{6,-2}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-6,62},{6,58}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Line(
          points={{0,78},{0,62}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{-94,78},{0,78}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{0,-80},{0,-94}},
          color={0,0,0},
          thickness=0.5),
        Line(
          points={{-94,50},{-12,50}},
          color={255,0,255},
          thickness=0.5),
        Line(
          points={{-94,-10},{-12,-10}},
          color={255,0,255},
          thickness=0.5),
        Line(
          points={{-94,-70},{-12,-70}},
          color={255,0,255},
          thickness=0.5),
        Line(
          points={{-16,54},{-12,50},{-16,46}},
          color={255,0,255},
          thickness=0.5),
        Line(
          points={{-16,-6},{-12,-10},{-16,-14}},
          color={255,0,255},
          thickness=0.5),
        Line(
          points={{-16,-66},{-12,-70},{-16,-74}},
          color={255,0,255},
          thickness=0.5),
        Text(
          extent={{-84,92},{80,82}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="Step Sequence")}),                         Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end InternalStepSequence_StateGraph1;
