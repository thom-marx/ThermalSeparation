within ThermalSeparation.Components.Columns;
package OldStartUpExample

  model MA_EQ_Cond_startUp "for different reflux ratios"
    import ThermalSeparation;
    //Lauftest: 10.12.2011

    inner ThermalSeparation.SystemTS systemTS(T_ref=298.15)
      annotation (Placement(transformation(extent={{60,100},{80,120}})));

    ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink1
      annotation (Placement(transformation(extent={{26,0},{46,20}}, rotation=0)));
    ThermalSeparation.Components.Columns.StructuredPackedColumn reactor(
      mapping={{1,1},{2,2},{3,3},{4,4}},
      redeclare package MediumVapour =
          ThermalSeparation.Media.Methylacetatsynthese_Vap,
      redeclare package MediumLiquid =
          ThermalSeparation.Media.Methylacetatsynthese_Liq,
      nS=4,
      stageLiquidFeed={5},
      x_l_start_in={0.1,0.3,0.2,0.4},
      x_l_start_out={0.7,0.1,0.1,0.1},
      x_v_start_in={0.1,0.3,0.2,0.4},
      x_v_start_out={0.7,0.1,0.1,0.1},
      x_l_profile=false,
      x_v_profile=false,
      T_l_profile=false,
      T_v_profile=false,
      hasLiquidFeed=false,
      wettedInitial=false,
      x_v_start_const={1e-5,1e-5,0.98,0.015},
      x_total_start={1e-1,1e-1,1 - 3e-1,1e-1},
      redeclare model Holdup =
          ThermalSeparation.Holdup.StructuredPackedColumn.Rocha,
      redeclare model HeatTransferWall = ThermalSeparation.Wall.Adiabatic (stat=
              true),
      redeclare model InitOption =
          ThermalSeparation.Components.Columns.BaseClasses.Initialization.Init_EquilibriumFilm,
      eps_liq_start=0.0045,
      considerStartUp=true,
      x_l_start_const={1e-5,1e-5,1 - 2e-5,0},
      k=0.001,
      redeclare model PressureLoss =
          ThermalSeparation.PressureLoss.StructuredPackedColumn.NominalLinear (
            Vdot_nom=0.0085, deltaP_nom=500),
      redeclare record Geometry =
          ThermalSeparation.Geometry.StructuredPackedColumn.MultipakI (
            hasInsulation=false),
      n_elements=8,
      redeclare model ThermoEquilibrium =
          ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid,
      redeclare model Reaction =
          ThermalSeparation.Reaction.ReactionKinetics.PseudoHomegeneousCatalysis
          (                                                                                              m_cat=1850),
      p_v_start_inlet=100000,
      p_v_start_outlet=100000,
      T_vap_start_bottom=373.15,
      T_vap_start_top=348.15,
      T_liq_start_bottom=373.15,
      T_liq_start_top=348.15,
      T_vapour_start=298.15,
      T_liquid_start=298.15,
      redeclare model BalanceEquations =
          ThermalSeparation.BalanceEquations.StructuredPackedColumn.NonEquilibrium.TwoPhaseVarState
          (                                                                                                                         redeclare
            model                                                                                                                                   FilmModel =
              ThermalSeparation.FilmModel.StructuredPackedColumn.MS))
                             annotation (Placement(transformation(extent={{-46,-18},{8,34}},
                            rotation=0)));

    ThermalSeparation.Components.SourcesSinks.SourceLiquid
      sourceLiquid_x(
      redeclare package MediumLiquid =
          ThermalSeparation.Media.Methylacetatsynthese_Liq,
      x={0,1,0,0},
      use_Flow=true,
      T=329.15)      "sinkGas.gasPortIn.Vdot/1000;sinkGas.gasPortIn.x"
      annotation (Placement(transformation(extent={{-66,66},{-46,86}})));
    ThermalSeparation.Components.SourcesSinks.CombLiquid_x2
      combLiquid_x_richtig1(
      redeclare package Medium =
          ThermalSeparation.Media.Methylacetatsynthese_Liq,
      x_start={0,0,1,0},
      V=1e-6,
      T_start=338.15)
      annotation (Placement(transformation(extent={{-22,64},{-2,84}})));
    Modelica.Blocks.Sources.Ramp Vdot_l(
      duration=30,
      height=3*1.43e-6,
      startTime=start_Vdot)
      annotation (Placement(transformation(extent={{-100,68},{-80,88}})));
    ThermalSeparation.Components.Columns.StructuredPackedColumn separator(
      k=0.001,
      mapping={{1,1},{2,2},{3,3},{4,4}},
      redeclare package MediumVapour =
          ThermalSeparation.Media.Methylacetatsynthese_Vap,
      redeclare package MediumLiquid =
          ThermalSeparation.Media.Methylacetatsynthese_Liq,
      nS=4,
      stageLiquidFeed={5},
      x_l_start_in={0.1,0.3,0.2,0.4},
      x_l_start_out={0.7,0.1,0.1,0.1},
      x_v_start_in={0.1,0.3,0.2,0.4},
      x_v_start_out={0.7,0.1,0.1,0.1},
      x_l_profile=false,
      x_v_profile=false,
      T_l_profile=false,
      T_v_profile=false,
      hasLiquidFeed=false,
      wettedInitial=false,
      redeclare model PressureLoss =
          ThermalSeparation.PressureLoss.StructuredPackedColumn.NominalLinear (
                                                              Vdot_nom=0.0085),
      redeclare record Geometry =
          ThermalSeparation.Geometry.StructuredPackedColumn.Rombopak_6M,
      x_v_start_const={1e-5,1e-5,1 - 3e-5,1e-5},
      x_total_start={3e-1,3e-1,1 - 3e-1*3,3e-1},
      x_l_start_const={1e-5,1e-5,1,1e-5},
      redeclare model Holdup =
          ThermalSeparation.Holdup.StructuredPackedColumn.Rocha,
      redeclare model HeatTransferWall = ThermalSeparation.Wall.Adiabatic (
            alpha_inner=20, stat=true),
      n_elements=4,
      redeclare model Reaction = ThermalSeparation.Reaction.NoReaction,
      redeclare model InitOption =
          ThermalSeparation.Components.Columns.BaseClasses.Initialization.Init_EquilibriumFilm,
      considerStartUp=true,
      eps_liq_start=(0.045 - 0.015),
      p_v_start_inlet=100000,
      p_v_start_outlet=100000,
      T_vap_start_bottom=373.15,
      T_vap_start_top=348.15,
      T_liq_start_bottom=373.15,
      T_liq_start_top=348.15,
      T_vapour_start=298.15,
      T_liquid_start=298.15,
      redeclare model BalanceEquations =
          ThermalSeparation.BalanceEquations.StructuredPackedColumn.NonEquilibrium.TwoPhaseVarState
          (redeclare model FilmModel =
              ThermalSeparation.FilmModel.StructuredPackedColumn.MS (redeclare
                replaceable model                                                                StateSelection =
                  ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.StateSelection1)))
                             annotation (Placement(transformation(extent={{-42,96},
              {12,148}},    rotation=0)));

    ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink2
      annotation (Placement(transformation(extent={{14,114},{34,134}},
                                                                    rotation=0)));
  output ThermalSeparation.Components.LiquidVolumes.ResultsLiquidVolumes result=tank.results;
  output ThermalSeparation.Components.Columns.BaseClasses.Results resultSump = sump.results;
  output ThermalSeparation.Components.Columns.BaseClasses.Results resultReactor = reactor.results;
  output ThermalSeparation.Components.Columns.BaseClasses.Results resultSeparator = separator.results;
  output Real PID_y = PID.y;
  output Real refluxRatio= refluxConverter.y;
  /*** monitoring ***/
  output Modelica.SIunits.MoleFraction x_1m[4]=reactor.x_l[2, :];
   Modelica.SIunits.MoleFraction x_2m[4]=reactor.x_l_in;
   Modelica.SIunits.MoleFraction x_3m[4]=separator.x_l_in;
   Modelica.SIunits.MoleFraction x_reboiler[4]=sump.x_l[1,:];
   Modelica.SIunits.Temperature T_1m=reactor.T_l[2];
   Modelica.SIunits.Temperature T_3m=separator.T_l_in;
   Modelica.SIunits.Temperature T_reboiler=sump.T[1];
    ThermalSeparation.Components.Condenser.TotalCondenser kondensator(
      redeclare package MediumVapour =
          ThermalSeparation.Media.Methylacetatsynthese_Vap,
      redeclare package MediumLiquid =
          ThermalSeparation.Media.Methylacetatsynthese_Liq,
      nS=4,
      mapping={{1,1},{2,2},{3,3},{4,4}},
      outletTempOption=ThermalSeparation.Components.Condenser.Enumerations.OutletTempOption.T_in,
      T_set=343.15)
      annotation (Placement(transformation(extent={{-17,-19},{17,19}},
          rotation=90,
          origin={43,187})));
    ThermalSeparation.Components.SourcesSinks.SplitLiquid_1p
      splitLiquid_x_2p(redeclare package Medium =
          ThermalSeparation.Media.Methylacetatsynthese_Liq, splitOption=
          ThermalSeparation.Components.SourcesSinks.Enumerations.SplitOption.out1_over_in)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
          origin={96,174})));
    ThermalSeparation.Components.LiquidVolumes.Tank         tank(
      redeclare package MediumLiquid =
          ThermalSeparation.Media.Methylacetatsynthese_Liq,
      nS=4,
      level_start=1e-5,
      x_l_start={1e-5,1e-5,1e-5,1 - 3e-5})
      annotation (Placement(transformation(extent={{76,36},{96,56}})));
    ThermalSeparation.Components.SourcesSinks.SinkLiquid          sinkLiquid1(
        redeclare package Medium =
          ThermalSeparation.Media.Methylacetatsynthese_Liq,
      use_Vdot=true,
      use_p=false,
      p=101000)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                                                                     rotation=270,
          origin={84,8})));
    Modelica.Blocks.Sources.Ramp q_flow2(
      startTime=1000,
      duration=30,
      height=0)
      annotation (Placement(transformation(extent={{54,-30},{74,-10}})));
    ThermalSeparation.Utilities.RefluxConverter refluxConverter
      annotation (Placement(transformation(extent={{74,86},{94,106}})));
    Modelica.Blocks.Sources.TimeTable switch(table=[0,10; 2400,10; 2450,1.7;
          800000,1.7])
      annotation (Placement(transformation(extent={{26,110},{46,130}})));
    ThermalSeparation.Components.Columns.StructuredPackedColumn sump(
      k=0.001,
      mapping={{1,1},{2,2},{3,3},{4,4}},
      redeclare package MediumVapour =
          ThermalSeparation.Media.Methylacetatsynthese_Vap,
      redeclare package MediumLiquid =
          ThermalSeparation.Media.Methylacetatsynthese_Liq,
      nS=4,
      stageLiquidFeed={5},
      x_l_start_in={0.1,0.3,0.2,0.4},
      x_l_start_out={0.7,0.1,0.1,0.1},
      x_v_start_in={0.1,0.3,0.2,0.4},
      x_v_start_out={0.7,0.1,0.1,0.1},
      x_l_profile=false,
      x_v_profile=false,
      T_l_profile=false,
      T_v_profile=false,
      hasLiquidFeed=false,
      wettedInitial=false,
      x_total_start={1e-1,1e-1,1 - 3e-1,1e-1},
      redeclare model PressureLoss =
          ThermalSeparation.PressureLoss.StructuredPackedColumn.NominalLinear (
            Vdot_nom=0.0085),
      redeclare model InitOption =
          ThermalSeparation.Components.Columns.BaseClasses.Initialization.Init_EquilibriumFilm,
      redeclare model Holdup =
          ThermalSeparation.Holdup.StructuredPackedColumn.NoOutflow,
      n_elements=1,
      redeclare model HeatTransferWall = ThermalSeparation.Wall.ConstQdot (
          alpha_inner=20,
          Qdot_set=6000,
          stat=true),
      x_v_start_const={1e-5,1e-5,1 - 3e-5,1e-5},
      considerStartUp=true,
      eps_liq_start=0.5*0.63,
      redeclare record Geometry =
          ThermalSeparation.Geometry.StructuredPackedColumn.MultipakI (
          hasInsulation=false,
          d=0.357,
          eps=0.99,
          H=2*0.83),
      redeclare model Reaction = ThermalSeparation.Reaction.NoReaction,
      x_l_start_const={1e-5,1e-5,1 - 2e-5,0},
      p_v_start_inlet=100000,
      p_v_start_outlet=100000,
      T_vap_start_bottom=373.15,
      T_vap_start_top=348.15,
      T_liq_start_bottom=373.15,
      T_liq_start_top=348.15,
      T_vapour_start=298.15,
      T_liquid_start=298.15,
      redeclare model BalanceEquations =
          ThermalSeparation.BalanceEquations.StructuredPackedColumn.NonEquilibrium.TwoPhaseVarState
          (redeclare model FilmModel =
              ThermalSeparation.FilmModel.StructuredPackedColumn.MS (redeclare
                replaceable model                                                                StateSelection =
                  ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.StateSelection3)))
                             annotation (Placement(transformation(extent={{-40,-88},
              {14,-36}},    rotation=0)));

    ThermalSeparation.Components.SourcesSinks.SourceGas          sourceGas_Vdot(
      redeclare package Medium =
          ThermalSeparation.Media.Methylacetatsynthese_Vap,
      x={1e-6,1e-6,1 - 3e-6,1e-6},
      use_Flow=false,
      Flow=0,
      T=338.65)
      annotation (Placement(transformation(extent={{-74,-96},{-54,-76}})));
    ThermalSeparation.Components.SourcesSinks.SinkLiquid          sinkLiquid(
        redeclare package Medium =
          ThermalSeparation.Media.Methylacetatsynthese_Liq)
      annotation (Placement(transformation(extent={{36,-94},{56,-74}})));
    ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink3
      annotation (Placement(transformation(extent={{32,-66},{52,-46}},
                                                                    rotation=0)));
    Modelica.Blocks.Sources.Ramp q_flow4(
      duration=30,
      height=0,
      startTime=0,
      offset=0.75)
      annotation (Placement(transformation(extent={{-18,44},{2,64}})));
    Modelica.Blocks.Continuous.LimPID PID(
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=1,
      initType=Modelica.Blocks.Types.InitPID.InitialOutput,
      y_start=10,
      k=1e-1,
      yMin=0.2,
      yMax=100) annotation (Placement(transformation(extent={{28,36},{48,56}})));
    Modelica.Blocks.Sources.TimeTable switch1(table=[0,1e6; 150000,1e6; 151000,2.5;
          1800e3,2.5])
      annotation (Placement(transformation(extent={{28,74},{48,94}})));

  parameter Real start_Reflux = 2400;
  parameter Real start_Vdot = 10;
    Modelica.Blocks.Sources.TimeTable Vdot(table=[0,0; start_Vdot,0; start_Vdot +
          30,5*1.43e-6; start_Vdot + 1000,5*1.43e-6; start_Vdot + 1030,1*1.43e-6;
          10000000,1*1.43e-6])
      annotation (Placement(transformation(extent={{-102,38},{-82,58}})));
  equation
       PID.u_m = if time < start_Reflux then PID.u_s else separator.x_l_in[1];
  //refluxConverter.u = if time < start_Reflux then switch1.y else PID.y;

    connect(separator.heatPort, ambientHeatSink2.heatPort1)     annotation (
        Line(
        points={{7.68,122},{8.87,122},{8.87,124},{12.4,124}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(q_flow2.y, sinkLiquid1.Vdot_set) annotation (Line(
        points={{75,-20},{70,-20},{70,-0.5},{92.5,-0.5}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(refluxConverter.y, splitLiquid_x_2p.split) annotation (Line(
        points={{96.1,96.9},{120,96.9},{120,198},{90,198},{90,174}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(q_flow4.y,PID. u_s) annotation (Line(
        points={{3,54},{14,54},{14,46},{26,46}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(switch.y, refluxConverter.u) annotation (Line(
        points={{47,120},{60.5,120},{60.5,97.2},{74,97.2}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(Vdot.y, sourceLiquid_x.Flow_In) annotation (Line(
        points={{-81,48},{-72,48},{-72,80},{-66,80}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(separator.upStreamOut, kondensator.vapourIn) annotation (Line(
        points={{-33.9,145.4},{-33.9,185.98},{26.28,185.98}},
        color={255,127,39},
        thickness=1,
        smooth=Smooth.None));
    connect(reactor.heatPort, ambientHeatSink1.heatPort1) annotation (Line(
        points={{3.68,8},{18.84,8},{18.84,10},{24.4,10}},
        color={188,51,69},
        smooth=Smooth.None));
    connect(reactor.downStreamOut, sump.downStreamIn) annotation (Line(
        points={{-0.1,-15.4},{-0.1,-28},{5.9,-28},{5.9,-38.6}},
        color={153,217,234},
        thickness=1,
        smooth=Smooth.None));
    connect(reactor.upStreamIn, sump.upStreamOut) annotation (Line(
        points={{-37.9,-15.4},{-37.9,-28},{-31.9,-28},{-31.9,-38.6}},
        color={255,127,39},
        thickness=1,
        smooth=Smooth.None));
    connect(ambientHeatSink3.heatPort1, sump.heatPort) annotation (Line(
        points={{30.4,-56},{20,-56},{20,-62},{9.68,-62}},
        color={188,51,69},
        smooth=Smooth.None));
    connect(sinkLiquid.liquidPortIn, sump.downStreamOut) annotation (Line(
        points={{36.4,-84},{22,-84},{22,-85.4},{5.9,-85.4}},
        color={153,217,234},
        thickness=1,
        smooth=Smooth.None));
    connect(sump.upStreamIn, sourceGas_Vdot.gasPortOut) annotation (Line(
        points={{-31.9,-85.4},{-41.95,-85.4},{-41.95,-86},{-52.6,-86}},
        color={255,127,39},
        thickness=1,
        smooth=Smooth.None));
    connect(combLiquid_x_richtig1.liquidPortIn2, sourceLiquid_x.liquidPortOut)
      annotation (Line(
        points={{-22,69},{-34,69},{-34,68},{-44.6,68},{-44.6,76}},
        color={153,217,234},
        thickness=1,
        smooth=Smooth.None));
    connect(separator.upStreamIn, reactor.upStreamOut) annotation (Line(
        points={{-33.9,98.6},{-33.9,31.4},{-37.9,31.4}},
        color={255,127,39},
        thickness=1,
        smooth=Smooth.None));
    connect(combLiquid_x_richtig1.liquidPortIn1, separator.downStreamOut)
      annotation (Line(
        points={{-22,79},{-22,88},{3.9,88},{3.9,98.6}},
        color={153,217,234},
        thickness=1,
        smooth=Smooth.None));
    connect(combLiquid_x_richtig1.liquidPortOut, reactor.downStreamIn)
      annotation (Line(
        points={{-2,74},{-0.1,74},{-0.1,31.4}},
        color={153,217,234},
        thickness=1,
        smooth=Smooth.None));
    connect(kondensator.liquidOut, splitLiquid_x_2p.liquidPortIn) annotation (
        Line(
        points={{59.72,185.98},{77.86,185.98},{77.86,184},{96,184}},
        color={153,217,234},
        thickness=1,
        smooth=Smooth.None));
    connect(tank.portIn, splitLiquid_x_2p.liquidPortOut2) annotation (Line(
        points={{86,55.6},{86,164},{91,164}},
        color={153,217,234},
        thickness=1,
        smooth=Smooth.None));
    connect(separator.downStreamIn, splitLiquid_x_2p.liquidPortOut1) annotation (
        Line(
        points={{3.9,145.4},{101,145.4},{101,164}},
        color={153,217,234},
        thickness=1,
        smooth=Smooth.None));
    connect(sinkLiquid1.liquidPortIn, tank.portOut) annotation (Line(
        points={{84,17.6},{86,17.6},{86,36.4}},
        color={153,217,234},
        thickness=1,
        smooth=Smooth.None));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              200}})),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
<p><ul>
<li> tray column</li>
<li> 3 components</li>
</ul></p>
</html>"),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,200}})));
  end MA_EQ_Cond_startUp;
end OldStartUpExample;
