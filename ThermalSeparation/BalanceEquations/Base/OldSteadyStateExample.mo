within ThermalSeparation.BalanceEquations.Base;
package OldSteadyStateExample

  model ProPen_RB_SteadyState
    "separation of alkanes, using a packed column, steady state balance equations"
    //Lauftest: 1.2.
    //Lauftest: 25.1.
    //Lauftest: 12.1. (3s, Kevin)

  ThermalSeparation.Components.SourcesSinks.SourceLiquid
                                                sourceLiquid(
      redeclare package MediumLiquid =
          ThermalSeparation.Media.Propane_Pentane_Liq,
      x={0.3,0.7},
      Flow=1,
      T=298.15)             annotation (Placement(transformation(extent={{-10,-10},
              {10,10}}, rotation=270,
          origin={8,88})));
  ThermalSeparation.Components.SourcesSinks.SinkLiquid
                                              sinkLiquid(redeclare package
        Medium =
          ThermalSeparation.Media.Propane_Pentane_Liq)
                        annotation (Placement(transformation(extent={{-10,-10},{
              10,10}},
                     rotation=270,
          origin={14,-32})));

  ThermalSeparation.Components.SourcesSinks.SinkGas
                                           sinkGas(redeclare package Medium =
          ThermalSeparation.Media.IdealGasMixtures.Propane_Pentane, p=300000)
                  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=90,
          origin={-30,84})));

    ThermalSeparation.Components.SourcesSinks.SourceGas
                                               sourceGas_Vdot(
      redeclare package Medium =
          ThermalSeparation.Media.IdealGasMixtures.Propane_Pentane,
      x={0.5,0.5},
      flowOption=ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowNdot,
      T=323.15)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
          origin={-28,-32})));

    Modelica.Blocks.Sources.Ramp ramp(
      duration=1,
      startTime=100,
      height=0,
      offset=450)
                annotation (Placement(transformation(extent={{-68,-64},{-48,-44}},
            rotation=0)));
    inner ThermalSeparation.SystemTS systemTS
      annotation (Placement(transformation(extent={{-74,-8},{-54,12}})));

  ThermalSeparation.Components.Columns.StructuredPackedColumn column1(
      hasLiquidFeed=false,
      hasVapourFeed=false,
      redeclare model Reaction = ThermalSeparation.Reaction.NoReaction,
      nS=2,
      mapping={{1,1},{2,2}},
      eps_liq_start=0.3,
      redeclare package MediumVapour =
          ThermalSeparation.Media.IdealGasMixtures.Propane_Pentane,
      redeclare package MediumLiquid =
          ThermalSeparation.Media.Propane_Pentane_Liq,
      redeclare record Geometry =
          ThermalSeparation.Geometry.StructuredPackedColumn.Geometry (H=10, d=4.9),
      redeclare model PressureLoss =
          ThermalSeparation.PressureLoss.StructuredPackedColumn.NominalLinear,
      redeclare model HeatTransferWall = ThermalSeparation.Wall.Adiabatic (stat=
              true),
      redeclare model ThermoEquilibrium =
          ThermalSeparation.PhaseEquilibrium.ConstantK (K_user={10,0.22}),
      redeclare model BalanceEquations =
          ThermalSeparation.BalanceEquations.StructuredPackedColumn.NonEquilibrium.TwoPhaseSteadyState,
      n_elements=1,
      wettedInitial=false,
      x_l_start_const={0.5,0.5},
      x_v_start_const={0.3,0.7},
      p_v_start_inlet=302000,
      p_v_start_outlet=300000,
      T_vapour_start=303.15,
      T_liquid_start=303.15,
      redeclare model InitOption =
          Components.Columns.BaseClasses.Initialization.None)
                           annotation (Placement(transformation(extent={{-34,18},
              {14,64}}, rotation=0)));

    ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink1
      annotation (Placement(transformation(extent={{36,32},{56,52}},rotation=0)));
  equation
    connect(ramp.y, sourceGas_Vdot.Flow_in) annotation (Line(
        points={{-47,-54},{-44,-54},{-44,-42},{-32,-42}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(ambientHeatSink1.heatPort1, column1.heatPort) annotation (Line(
        points={{34.4,42},{22,42},{22,41},{10.16,41}},
        color={188,51,69},
        smooth=Smooth.None));
    connect(sourceGas_Vdot.gasPortOut, column1.upStreamIn) annotation (Line(
        points={{-28,-20.6},{-28,20.3},{-26.8,20.3}},
        color={255,127,39},
        thickness=1,
        smooth=Smooth.None));
    connect(sinkLiquid.liquidPortIn, column1.downStreamOut) annotation (Line(
        points={{14,-22.4},{12,-22.4},{12,20.3},{6.8,20.3}},
        color={153,217,234},
        thickness=1,
        smooth=Smooth.None));
    connect(column1.upStreamOut, sinkGas.gasPortIn) annotation (Line(
        points={{-26.8,61.7},{-26.8,67.85},{-30,67.85},{-30,74.4}},
        color={255,127,39},
        thickness=1,
        smooth=Smooth.None));
    connect(column1.downStreamIn, sourceLiquid.liquidPortOut) annotation (Line(
        points={{6.8,61.7},{6.8,67.85},{8,67.85},{8,76.6}},
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
  end ProPen_RB_SteadyState;
end OldSteadyStateExample;
