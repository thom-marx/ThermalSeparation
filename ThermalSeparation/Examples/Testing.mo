within ThermalSeparation.Examples;
package Testing
  package Complex
    model CycleComplex
      //import ThermalSeparation;
      Real x_CO2_removed = max(1e-7, 1 - Absorber.Ndot_v[Absorber.n] * Absorber.x_v[Absorber.n, 3] / (Absorber.Ndot_v_in * Absorber.x_v_in[3]));
      Modelica.Units.SI.MassFlowRate water_diff = Absorber.Ndot_v_in * Absorber.x_v_in[1] * 0.018 - Absorber.Ndot_v[Absorber.n] * Absorber.x_v[Absorber.n, 1] * 0.018;
      //Real qspec = sumpV.Q_in/(1000^2) * 1/(max(0.04,Desorber.mdot_v[Desorber.n]*Desorber.X_v[Desorber.n,2]));
      // Modelica.Units.SI.MassFlowRate water_diff1_gas = Desorber.Ndot_v_in*Desorber.x_v_in[1]*0.018 - Desorber.Ndot_v[Desorber.n]*Desorber.x_v[Desorber.n,1]*0.018;
      // Modelica.Units.SI.MassFlowRate water_diff1_liq = Desorber.Ndot_l_in*Desorber.x_l_in[1]*0.018 - Desorber.Ndot_l[1]*Desorber.x_l[1,1]*0.018;
      //
      // Modelica.Units.SI.MassFlowRate bal_reboiler = sumpV.liquidPortIn.Vdot*sumpV.mediumLiquidIn.d + sumpV.gasPortOut.Vdot*sumpV.mediumVapour.d + sumpV.liquidPortOut.Vdot*sumpV.mediumLiquid.d;
      // Modelica.Units.SI.MassFlowRate bal_condenser = flash_Condenser_simple.gasPortIn.Vdot*flash_Condenser_simple.mediumVapourIn.d + flash_Condenser_simple.gasPortOut.Vdot*flash_Condenser_simple.mediumVapour.d + flash_Condenser_simple.liquidPortOut.Vdot*flash_Condenser_simple.mediumLiquid.d;
      // Modelica.Units.SI.MassFlowRate bal_desorber = Desorber.upStreamIn.Vdot*Desorber.mediumVapourIn.d + Desorber.downStreamIn.Vdot*Desorber.mediumLiquidIn.d+Desorber.upStreamOut.Vdot*Desorber.mediumVapour[20].d+Desorber.downStreamOut.Vdot*Desorber.mediumLiquid[1].d;
      // Modelica.Units.SI.MassFlowRate bal_absorber = Absorber.upStreamIn.Vdot*Absorber.mediumVapourIn.d + Absorber.downStreamIn.Vdot*Absorber.mediumLiquidIn.d+Absorber.upStreamOut.Vdot*Absorber.mediumVapour[30].d+Absorber.downStreamOut.Vdot*Absorber.mediumLiquid[1].d;
      // Modelica.Units.SI.MassFlowRate bal_pump = idealPumpControlledVdot.liquidIn.Vdot*idealPumpControlledVdot.mediumLiquidIn.d+idealPumpControlledVdot.liquidOut.Vdot*idealPumpControlledVdot.mediumLiquid.d;
      // Modelica.Units.SI.MassFlowRate bal_pump1 = idealPumpControlledVdot1.liquidIn.Vdot*idealPumpControlledVdot1.mediumLiquidIn.d+idealPumpControlledVdot1.liquidOut.Vdot*idealPumpControlledVdot1.mediumLiquid.d;
      // Modelica.Units.SI.MassFlowRate bal_HEX = counterFlowHeatExchanger.hotLiquidIn.Vdot*counterFlowHeatExchanger.hotLiquidTube.mediumLiquidIn.d+counterFlowHeatExchanger.coldLiquidIn.Vdot*counterFlowHeatExchanger.coldLiquidTube.mediumLiquidIn.d+counterFlowHeatExchanger.hotLiquidOut.Vdot*counterFlowHeatExchanger.hotLiquidTube.mediumLiquid.d+counterFlowHeatExchanger.coldLiquidOut.Vdot*counterFlowHeatExchanger.coldLiquidTube.mediumLiquid.d;
      // Modelica.Units.SI.MassFlowRate bal_cooler = cooler.hotLiquidIn.Vdot*cooler.mediumLiquidIn.d + cooler.coldLiquidOut.Vdot*cooler.mediumLiquid.d;
      // Modelica.Units.SI.MassFlowRate bal_tank = tank.portIn.Vdot*tank.mediumLiquidIn.d+tank.portOut.Vdot*tank.mediumLiquid.d;
      // Modelica.Units.SI.MassFlowRate bal_tank1 = tank1.portIn.Vdot*tank1.mediumLiquidIn.d+tank1.portOut.Vdot*tank1.mediumLiquid.d;
      // Modelica.Units.SI.MassFlowRate bal_comb = combLiquid_x.liquidPortIn1.Vdot*combLiquid_x.mediumIn1.d+combLiquid_x.liquidPortIn2.Vdot*combLiquid_x.mediumIn2.d+combLiquid_x.liquidPortOut.Vdot*combLiquid_x.mediumOut.d;
      // Modelica.Units.SI.MassFlowRate bal_comb1 = combLiquid_x1.liquidPortIn1.Vdot*combLiquid_x1.mediumIn1.d+combLiquid_x1.liquidPortIn2.Vdot*combLiquid_x1.mediumIn2.d+combLiquid_x1.liquidPortOut.Vdot*combLiquid_x1.mediumOut.d;
      //
      // Modelica.Units.SI.MassFlowRate bal_total = Absorber.upStreamIn.Vdot*Absorber.mediumVapourIn.d+Absorber.upStreamOut.Vdot*Absorber.mediumVapour[30].d+ flash_Condenser_simple.gasPortOut.Vdot*flash_Condenser_simple.mediumVapour.d+combLiquid_x1.liquidPortIn1.Vdot*combLiquid_x1.mediumIn1.d;
      ThermalSeparation.Components.Columns.StructuredPackedColumn Absorber(inertVapour = {false, true, false, true}, inertLiquid = {false, false, true}, nS = 2, redeclare model HeatTransferWall =
            ThermalSeparation.Wall.Adiabatic,                                                                                                                                                                                         mapping = {{1, 1}, {3, 2}}, wettedInitial = false, x_l_start_const = {0.88, 0.01, 0.11}, redeclare record Geometry =
            ThermalSeparation.Geometry.StructuredPackedColumn.Mellapak250Y (                                                                                                                                                                                                        d = 15, H = 15), redeclare model PressureLoss =
            ThermalSeparation.PressureLoss.StructuredPackedColumn.Stichlmair,                                                                                                                                                                                                        redeclare model Holdup =
            ThermalSeparation.Holdup.StructuredPackedColumn.StichlmairStat,                                                                                                                                                                                                        eps_liq_start = 0.05, x_v_start_const = {0.1841, 0.0471, 0.08, 0.7}, redeclare package MediumVapour =
            Media.H2O_O2_CO2_N2_Vap,                                                                                                                                                                                                        redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare model ThermoEquilibrium =
            PhaseEquilibrium.H2O_CO2_MEA (                                                                                                                                                                                                        factor_K = {1, 1.1, 1}), redeclare model Reaction =
            Reaction.ReactionEquilibrium.CO2_MEA,
        redeclare model BalanceEquations = BalanceEquations.StructuredPackedColumn.NonEquilibrium.TwoPhaseVarState (redeclare model FilmModel = ThermalSeparation.FilmModel.StructuredPackedColumn.TrueEquilibrium (redeclare model StateSelection = ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.StateSelection3)),
                                                                                                                                                                                                        redeclare model InitOption =
            Components.Columns.BaseClasses.Initialization.Dyncap1_GG,
        p_v_start_inlet=100500,
        p_v_start_outlet=100000,
        T_vapour_start=323.15,
        T_liquid_start=323.15,                                                                                                                                                                                                        n_elements = 25) annotation (
        Placement(transformation(extent = {{-88, -42}, {-42, 4}})));
      ThermalSeparation.Components.SourcesSinks.SinkGas sinkGas(redeclare package Medium =
            Media.H2O_O2_CO2_N2_Vap,                                                                                p = 97000) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-81, 46})));
      inner ThermalSeparation.SystemTS systemTS annotation (
        Placement(transformation(extent = {{-100, 80}, {-80, 100}})));
      ThermalSeparation.Components.SourcesSinks.SourceGas sourceGas_Vdot(flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowMdot, x = {0.075, 0.034, 0.145, 0.746}, redeclare package Medium =
            Media.H2O_O2_CO2_N2_Vap,                                                                                                                                                                                                        T = 308.15) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-81, -66})));
      Modelica.Blocks.Sources.Ramp Vdot(offset = 330, duration = 0.1, height = -0.2 * 330, startTime = 2000) annotation (
        Placement(transformation(extent = {{-77, -94}, {-85, -86}})));
      ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink annotation (
        Placement(transformation(extent = {{-38, -24}, {-26, -12}})));
      ThermalSeparation.Components.HeatExchanger.CcounterFlowHeatExchanger counterFlowHeatExchanger(c_l_start_hot = fill(1000, 3), c_l_start_cold = fill(1000, 3), ss_hot_c2 = false, ss_cold_c2 = false, nScold = 2, surfaceArea = 3, redeclare package MediumLiquidCold =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare package MediumLiquidHot =
            Media.H2O_CO2_MEA_Liq)                                                                                                                                                                                                         annotation (
        Placement(transformation(extent = {{-14, -2}, {20, -32}})));
      ThermalSeparation.Components.HeatExchanger.LiquidCooler cooler(T_set = 313.15, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq)                                                                                                           annotation (
        Placement(transformation(extent = {{7, -9}, {-7, 9}}, rotation = 270, origin = {-33, 15})));
      ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink1 annotation (
        Placement(transformation(extent = {{-20, 28}, {-12, 36}})));
      ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink2 annotation (
        Placement(transformation(extent = {{130, 15}, {142, 24}})));
      ThermalSeparation.Components.SourcesSinks.SinkGas sinkGas1(redeclare package Medium =
            Media.H2O_CO2_Vap,                                                                                 p = 200000) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {97, 86})));
      ThermalSeparation.Components.Reboiler.KettleReboilerEq sumpV(A = 50, H = 2, inert_Liquid = {false, false, true}, A_HT = 150, eps_liq_init = 0.58, x_total_start = {0.8, 0.2}, init_standalone = false, fixed_mol_init = 7, redeclare model InnerHT =
            ThermalSeparation.HeatAndMassTransfer.HTResistance.NoHTResistance,                                                                                                                                                                                                        redeclare package MediumVapour =
            Media.H2O_CO2_Vap,                                                                                                                                                                                                        redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare model ThermoEquilibrium =
            PhaseEquilibrium.H2O_CO2_MEA,                                                                                                                                                                                                        redeclare model Reaction =
            Reaction.ReactionEquilibrium.CO2_MEA,                                                                                                                                                                                                        init_option = ThermalSeparation.Components.Reboiler.InitOptionEq.init_mol, T_init = 397.15, p_init = 213000) annotation (
        Placement(transformation(extent = {{94, -38}, {114, -18}})));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation (
        Placement(transformation(extent = {{70, -35}, {84, -21}})));
      ThermalSeparation.Components.Columns.StructuredPackedColumn Desorber(mapping = {{1, 1}, {2, 2}}, inertLiquid = {false, false, true}, redeclare model HeatTransferWall =
            ThermalSeparation.Wall.Adiabatic,                                                                                                                                                                   T_l_profile = false, x_l_start_in = {0.8, 0.1, 0.1}, x_l_start_out = {0.5, 0.2, 0.3}, x_l_profile = false, wettedInitial = true, T_v_profile = false, x_l_start_const = {0.6, 0.25, 0.15}, x_v_start_const = {0.35, 0.65}, nS = 2, redeclare model PressureLoss =
            ThermalSeparation.PressureLoss.StructuredPackedColumn.Stichlmair,                                                                                                                                                                                                        redeclare model Holdup =
            ThermalSeparation.Holdup.StructuredPackedColumn.StichlmairStat,                                                                                                                                                                                                        eps_liq_start = 0.05, redeclare package MediumVapour =
            Media.H2O_CO2_Vap,                                                                                                                                                                                                        redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare model ThermoEquilibrium =
            PhaseEquilibrium.H2O_CO2_MEA,                                                                                                                                                                                                        redeclare model Reaction =
            Reaction.ReactionEquilibrium.CO2_MEA,                                                                                                                                                                                                        redeclare record Geometry =
            Geometry.StructuredPackedColumn.Mellapak250Y (                                                                                                                                                                                                        d = 8.65),
        redeclare model BalanceEquations = BalanceEquations.StructuredPackedColumn.NonEquilibrium.TwoPhaseVarState (redeclare model FilmModel = ThermalSeparation.FilmModel.StructuredPackedColumn.TrueEquilibrium (redeclare model StateSelection = ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.StateSelection3)),
                                                                                                                                                                                                        redeclare model InitOption =
            Components.Columns.BaseClasses.Initialization.Dyncap1_GG,                                                                                                                                                                                                        n_elements = 30,
        p_v_start_inlet=200500,
        p_v_start_outlet=200000,
        T_vap_start_bottom=383.15,
        T_vap_start_top=373.15,
        T_vapour_start=373.15,
        T_liquid_start=393.15)                                                                                                                                                                                                         annotation (
        Placement(transformation(extent = {{80, -5}, {126, 41}})));
      ThermalSeparation.Components.Pumps.idealPumpControlledVdot idealPumpControlledVdot(V_flow_start = 0.9, T = 1, dT = 0.5, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq)                                                                                                                                                    annotation (
        Placement(transformation(extent = {{11, -10}, {-11, 10}}, rotation = 90, origin = {51, -88})));
      Modelica.Blocks.Continuous.LimPID PID(yMin = 0, y_start = 0.9, yMax = 1.2, controllerType = Modelica.Blocks.Types.SimpleController.PI, initType = Modelica.Blocks.Types.Init.InitialOutput, k = 1, Ti = 10) annotation (
        Placement(transformation(extent = {{4, -4}, {-4, 4}}, rotation = 0, origin = {34, -54})));
      ThermalSeparation.Components.LiquidVolumes.Tank tank(inertLiquid = {false, false, true}, x_l_start = {0.868, 0.022, 0.11}, d_volume = 3, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                     p_gas = 100000, T_start = 393.15, level_start = 10) annotation (
        Placement(transformation(extent = {{8.5, -9}, {-8.5, 9}}, rotation = 0, origin = {102.5, -79})));
      ThermalSeparation.Components.LiquidVolumes.Tank tank1(inertLiquid = {false, false, true}, d_volume = 1, x_l_start = {0.84, 0.053, 0.107}, level_start = 10, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                        p_gas = 100000, T_start = 313.15) annotation (
        Placement(transformation(extent = {{-63, -68}, {-43, -48}})));
      ThermalSeparation.Components.Pumps.idealPumpControlledVdot idealPumpControlledVdot1(V_flow_start = 0.9, T = 1, dT = 0.5, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq)                                                                                                                                                     annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {-24, -88})));
      Modelica.Blocks.Continuous.LimPID PID1(controllerType = Modelica.Blocks.Types.SimpleController.PI, k = 100, Ti = 0.1, yMin = 0, initType = Modelica.Blocks.Types.Init.InitialOutput, y_start = 0.1, yMax = 1.2) annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 0, origin = {-26, -55})));
      Modelica.Blocks.Sources.RealExpression captureRatio1(y = x_CO2_removed) annotation (
        Placement(transformation(extent = {{-5, -4}, {5, 4}}, rotation = 180, origin = {43, -62})));
      Modelica.Blocks.Sources.RealExpression level(y = 10) annotation (
        Placement(transformation(extent={{-38,-72},{-28,-64}})));
      ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink3 annotation (
        Placement(transformation(extent = {{132, 56}, {144, 68}})));
      ThermalSeparation.Components.SourcesSinks.CombLiquid_x combLiquid_x(x_start = {0.868, 0.022, 0.11}, redeclare package Medium =
            Media.H2O_CO2_MEA_Liq,                                                                                                                          T_start = 393.15) annotation (
        Placement(transformation(extent = {{-5, -5}, {5, 5}}, rotation = 270, origin = {100, -45})));
      Modelica.Blocks.Continuous.LimPID PID2(yMax = 300000000, yMin = 50000000, y_start = 201000000, initType = Modelica.Blocks.Types.Init.InitialOutput, controllerType = Modelica.Blocks.Types.SimpleController.PI, k = 10000, Ti = 0.005) annotation (
        Placement(transformation(extent = {{-4, 4}, {4, -4}}, rotation = 0, origin = {62, -28})));
      Modelica.Blocks.Sources.RealExpression T_Reboiler(y = sumpV.T) annotation (
        Placement(transformation(extent = {{-6, -4}, {6, 4}}, rotation = 0, origin = {46, -16})));
      Modelica.Blocks.Sources.RealExpression T_set(y = 124 + 273) annotation (
        Placement(transformation(extent = {{40, -31}, {50, -25}})));
      ThermalSeparation.Components.SourcesSinks.CombLiquid_x combLiquid_x1(x_start = {0.868, 0.022, 0.11}, redeclare package Medium =
            Media.H2O_CO2_MEA_Liq,                                                                                                                           T_start = 393.15) annotation (
        Placement(transformation(extent = {{-5, -5}, {5, 5}}, rotation = 270, origin = {102, -59})));
      ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid_x(flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowMdot, Flow = 1021, use_Flow = true, x = {1, 0, 0}, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        T = 393.15) annotation (
        Placement(transformation(extent = {{138, -46}, {118, -66}})));
      Modelica.Blocks.Sources.RealExpression Water_in(y = -water_diff) annotation (
        Placement(transformation(extent = {{5, -4}, {-5, 4}}, rotation = 180, origin = {124, -82})));
      Modelica.Blocks.Sources.Ramp Vdot1(duration = 0.1, height = 0.05, offset = 0.85, startTime = 200) annotation (
        Placement(transformation(extent = {{66, -57}, {60, -51}})));
      Components.Condenser.FlashCondenser_CO2_H2O flashCondenser_CO2_H2O annotation (
        Placement(transformation(extent = {{90, 48}, {110, 68}})));
      Modelica.Blocks.Sources.Ramp Vdot2(duration = 0.1, height = 0.05, offset = 0.85, startTime = 200) annotation (
        Placement(transformation(extent = {{78, -77}, {72, -71}})));
    equation
      connect(Vdot.y, sourceGas_Vdot.Flow_in) annotation (
        Line(points = {{-85.4, -90}, {-96, -90}, {-96, -76}, {-85, -76}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(prescribedHeatFlow.port, sumpV.heatPort) annotation (
        Line(points = {{84, -28}, {112.2, -28}}, color = {191, 0, 0}, smooth = Smooth.None));
      connect(PID1.y, idealPumpControlledVdot1.VdotInput) annotation (
        Line(points = {{-21.6, -55}, {-20, -55}, {-20, -76}, {-24, -76}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(PID1.u_s, tank1.y) annotation (
        Line(points = {{-30.8, -55}, {-44, -55}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(level.y, PID1.u_m) annotation (
        Line(points={{-27.5,-68},{-26,-68},{-26,-59.8}},        color = {0, 0, 127}, smooth = Smooth.None));
      connect(T_Reboiler.y, PID2.u_m) annotation (
        Line(points = {{52.6, -16}, {62, -16}, {62, -23.2}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(T_set.y, PID2.u_s) annotation (
        Line(points = {{50.5, -28}, {57.2, -28}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(Water_in.y, sourceLiquid_x.Flow_In) annotation (
        Line(points = {{129.5, -82}, {138, -82}, {138, -60}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(PID2.y, prescribedHeatFlow.Q_flow) annotation (
        Line(points = {{66.4, -28}, {70, -28}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(PID.u_s, Vdot1.y) annotation (
        Line(points = {{38.8, -54}, {59.7, -54}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(captureRatio1.y, PID.u_m) annotation (
        Line(points = {{37.5, -62}, {34, -62}, {34, -58.8}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(cooler.heatPort, ambientHeatSink1.heatPort1) annotation (
        Line(points = {{-33, 20.74}, {-33, 32}, {-20.64, 32}}, color = {188, 51, 69}, thickness = 1, smooth = Smooth.None));
      connect(counterFlowHeatExchanger.hotLiquidOut, cooler.hotLiquidIn) annotation (
        Line(points = {{-14, -12.2}, {-14, 8}, {-29.94, 8}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(Absorber.downStreamIn, cooler.coldLiquidOut) annotation (
        Line(points = {{-48.9, 1.7}, {-48.9, 8}, {-35.88, 8}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(sinkGas.gasPortIn, Absorber.upStreamOut) annotation (
        Line(points = {{-81, 36.4}, {-81, 21.2}, {-81.1, 21.2}, {-81.1, 1.7}}, color = {255, 127, 39}, thickness = 1, smooth = Smooth.None));
      connect(sourceGas_Vdot.gasPortOut, Absorber.upStreamIn) annotation (
        Line(points = {{-81, -54.6}, {-81, -46.3}, {-81.1, -46.3}, {-81.1, -39.7}}, color = {255, 127, 39}, thickness = 1, smooth = Smooth.None));
      connect(tank1.portIn, Absorber.downStreamOut) annotation (
        Line(points = {{-53, -48.4}, {-53, -44.1}, {-48.9, -44.1}, {-48.9, -39.7}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(tank1.portOut, idealPumpControlledVdot1.liquidIn) annotation (
        Line(points = {{-53, -67.6}, {-53, -88}, {-33.4, -88}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(idealPumpControlledVdot1.liquidOut, counterFlowHeatExchanger.coldLiquidIn) annotation (
        Line(points = {{-13.8, -88}, {-10, -88}, {-10, -44}, {-14, -44}, {-14, -22.1}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(counterFlowHeatExchanger.hotLiquidIn, idealPumpControlledVdot.liquidOut) annotation (
        Line(points={{20,-22.4},{20,-22},{38,-22},{38,-88},{40.8,-88}},
                                                             color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(counterFlowHeatExchanger.coldLiquidOut, Desorber.downStreamIn) annotation (
        Line(points={{20,-12.2},{36,-12.2},{36,96},{119.1,96},{119.1,38.7}},            color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(combLiquid_x1.liquidPortIn2, combLiquid_x.liquidPortOut) annotation (
        Line(points = {{99.5, -54}, {100, -54}, {100, -50}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(sourceLiquid_x.liquidPortOut, combLiquid_x1.liquidPortIn1) annotation (
        Line(points = {{116.6, -56}, {116.6, -52}, {104.5, -52}, {104.5, -54}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(Desorber.heatPort, ambientHeatSink2.heatPort1) annotation (
        Line(points = {{122.32, 18}, {126, 18}, {126, 19.5}, {129.04, 19.5}}, color = {188, 51, 69}, thickness = 1, smooth = Smooth.None));
      connect(Absorber.heatPort, ambientHeatSink.heatPort1) annotation (
        Line(points = {{-45.68, -19}, {-43.84, -19}, {-43.84, -18}, {-38.96, -18}}, color = {188, 51, 69}, thickness = 1, smooth = Smooth.None));
      connect(sinkGas1.gasPortIn, flashCondenser_CO2_H2O.gasPortOut) annotation (
        Line(points = {{97, 76.4}, {97, 72.2}, {99.6, 72.2}, {99.6, 66.6}}, color = {255, 127, 39}, thickness = 1));
      connect(flashCondenser_CO2_H2O.gasPortIn, Desorber.upStreamOut) annotation (
        Line(points = {{90.8, 58}, {86.9, 58}, {86.9, 38.7}}, color = {255, 127, 39}, thickness = 1));
      connect(flashCondenser_CO2_H2O.heatPort, ambientHeatSink3.heatPort1) annotation (
        Line(points = {{108, 58}, {119, 58}, {119, 62}, {131.04, 62}}, color = {188, 51, 69}));
      connect(flashCondenser_CO2_H2O.liquidPortOut, combLiquid_x.liquidPortIn1) annotation (
        Line(points = {{99.5, 49.4}, {144, 49.4}, {144, -40}, {102.5, -40}}, color = {153, 217, 234}, thickness = 1));
      connect(combLiquid_x.liquidPortIn2, sumpV.liquidPortOut) annotation (
        Line(points = {{97.5, -40}, {100, -40}, {100, -36.8}, {103.2, -36.8}}, color = {153, 217, 234}, thickness = 1));
      connect(Desorber.downStreamOut, sumpV.liquidPortIn) annotation (
        Line(points = {{119.1, -2.7}, {119.1, -10.35}, {106.8, -10.35}, {106.8, -19.8}}, color = {153, 217, 234}, thickness = 1));
      connect(Desorber.upStreamIn, sumpV.gasPortOut) annotation (
        Line(points = {{86.9, -2.7}, {86.9, -10.35}, {100, -10.35}, {100, -19.8}}, color = {255, 127, 39}, thickness = 1));
      connect(combLiquid_x1.liquidPortOut, tank.portIn) annotation (
        Line(points = {{102, -64}, {102, -70.36}, {102.5, -70.36}}, color = {153, 217, 234}, thickness = 1));
      connect(tank.portOut, idealPumpControlledVdot.liquidIn) annotation (
        Line(points = {{102.5, -87.64}, {102.5, -94}, {72, -94}, {72, -88}, {60.4, -88}}, color = {153, 217, 234}, thickness = 1));
      connect(PID.y, idealPumpControlledVdot.VdotInput) annotation (
        Line(points = {{29.6, -54}, {26, -54}, {26, -56}, {26, -72}, {51, -72}, {51, -74.8}}, color = {0, 0, 127}));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {150, 100}})),
        experiment(StopTime=10000, __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput,
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {150, 100}})));
    end CycleComplex;

    model Cycle_Heyden_StartUp
      Real x_CO2_removed = max(1e-7, 1 - Absorber.x_v[Absorber.n, 3] / max(0.001, sourceGas_Vdot.x[3]));
      Modelica.Units.SI.MassFlowRate water_diff = Absorber.Ndot_v_in * Absorber.x_v_in[1] * 0.018 - Absorber.Ndot_v[Absorber.n] * Absorber.x_v[Absorber.n, 1] * 0.018;
      Real qspec = sumpV.Q_in / 1000 ^ 2 * 1 / max(0.04, Desorber.mdot_v[Desorber.n] * Desorber.X_v[Desorber.n, 2]);
      ThermalSeparation.Components.Columns.StructuredPackedColumn_newStartUpShutDown Absorber(inertVapour = {false, true, false, true}, inertLiquid = {false, false, true}, nS = 2, mapping = {{1, 1}, {3, 2}}, wettedInitial = false, redeclare model Holdup =
            ThermalSeparation.Holdup.StructuredPackedColumn.StichlmairStat,                                                                                                                                                                                                        considerStartUp = true, StartUp_CCS = true, friggelfaktor = 0.2e5, smooth_startUp = false, delay_startUp = 750, lowBoilingPoint = {false, true, true, true}, eps_liq_start = 0.001, x_v_start_const = {0.048, 0.2, 0.001, 0.751}, switchingCondition_Boiling = false, x_v_dummy = {1, 0, 0, 0}, y_PID = 1000, gain = 100, Vdot_startUp_pressure = 5, x_l_start_const = {0.8593, 0.0325, 0.1082}, redeclare package MediumVapour =
            ThermalSeparation.Media.H2O_O2_CO2_N2_Vap,                                                                                                                                                                                                        redeclare package MediumLiquid =
            ThermalSeparation.Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare model Reaction =
            Reaction.ReactionEquilibrium.CO2_MEA,                                                                                                                                                                                                        n_elements = 15, redeclare model FilmModel =
            FilmModel.StructuredPackedColumn.TrueEquilibriumStartUpCCSAbsorption (                                                                                                                                   faktor_Ndot_v = 30000, redeclare model StateSelection =
                ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.StateSelection3),                                                                                                                                                                                                        redeclare model PressureLoss =
            ThermalSeparation.PressureLoss.StructuredPackedColumn.NominalLinear (                                                                                                                                    deltaP_nom(displayUnit = "Pa") = 10, Vdot_nom = 1), redeclare model InitOption =
            Components.Columns.BaseClasses.Initialization.DyncapStartUpAbsorption,                                                                                                                                                                                                        redeclare record Geometry =
            ThermalSeparation.Geometry.StructuredPackedColumn.Mellapak250Y (                                                                                                                                                                                                        d = 14.56, H = 15), redeclare model HeatTransferWall =
            Wall.ConstAlpha,                                                                                                                                                                                                        redeclare model ThermoEquilibrium =
            PhaseEquilibrium.H2O_CO2_MEA_startUp (                                                                                                                                                                                                        factor_K = {1, 1.03, 1}), p_v_start_inlet = 100000, p_v_start_outlet = 100000, T_vapour_start = 308.15, T_liquid_start = 308.15) annotation (
        Placement(transformation(extent = {{-178, -28}, {-80, 70}})));
      ThermalSeparation.Components.SourcesSinks.SinkGas sinkGas(redeclare package Medium =
            Media.H2O_O2_CO2_N2_Vap,                                                                                p = 102600) annotation (
        Placement(transformation(extent = {{-17, -21.5}, {17, 21.5}}, rotation = 90, origin = {-162.5, 157})));
      ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink annotation (
        Placement(transformation(extent = {{-72, 14}, {-60, 26}})));
      ThermalSeparation.Components.HeatExchanger.CcounterFlowHeatExchanger counterFlowHeatExchanger(c_l_start_hot = fill(1000, 3), c_l_start_cold = fill(1000, 3), nScold = 2, ss_hot_c2 = true, ss_cold_c2 = true, alpha_input = false, surfaceArea = 550, mFlowHotLiquid_nom = 40, mFlowColdLiquid_nom = 40, redeclare package MediumLiquidCold =
            ThermalSeparation.Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare package MediumLiquidHot =
            ThermalSeparation.Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        ThotLiquidIn_start = 313.15, TcoldLiquidIn_start = 298.15, ThotLiquid_start = 313.15, TcoldLiquid_start = 298.15) annotation (
        Placement(transformation(extent = {{-24, 24}, {36, -30}})));
      ThermalSeparation.Components.HeatExchanger.LiquidCooler cooler(T_set = 313.15, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq)                                                                                                           annotation (
        Placement(transformation(extent = {{7, -9}, {-7, 9}}, rotation = 270, origin = {-77, 93})));
      ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink1 annotation (
        Placement(transformation(extent = {{-4, -4}, {4, 4}}, rotation = 90, origin = {-76, 110})));
      ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink2 annotation (
        Placement(transformation(extent = {{130, 15}, {142, 24}})));
      ThermalSeparation.Components.SourcesSinks.SinkGas sinkGas1(redeclare package Medium =
            Media.H2O_CO2_Vap,                                                                                 p = 200000) annotation (
        Placement(transformation(extent = {{-16, -17.5}, {16, 17.5}}, rotation = 90, origin = {146.5, 164})));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation (
        Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 180, origin = {167, -102})));
      ThermalSeparation.Components.Columns.StructuredPackedColumn_newStartUpShutDown Desorber(mapping = {{1, 1}, {2, 2}}, inertLiquid = {false, false, true}, T_l_profile = false, x_l_start_in = {0.8, 0.1, 0.1}, x_l_start_out = {0.5, 0.2, 0.3}, x_l_profile = false, T_v_profile = false, nS = 2, inertVapour = {false, false}, wettedInitial = false, eps_liq_start = 0.001, considerStartUp = true, StartUp_CCS = true, switchingCondition_Boiling = true, lowBoilingPoint = {false, true}, redeclare model InitOption =
            ThermalSeparation.Components.Columns.BaseClasses.Initialization.DyncapStartUpDesorption,                                                                                                                                                                                                        redeclare model Holdup =
            ThermalSeparation.Holdup.StructuredPackedColumn.StichlmairStat,                                                                                                                                                                                                        x_v_start_const = {0, 1}, x_l_start_const = {0.8593, 0.0325, 0.1082}, y_PID = 1000, Vdot_startUp_pressure = 5, gain = 100, smooth_startUp = true, k = 0.04, delay_startUp = 100, friggelfaktor = 0.7e5, redeclare model Reaction =
            Reaction.ReactionEquilibrium.CO2_MEA,                                                                                                                                                                                                        redeclare record Geometry =
            Geometry.StructuredPackedColumn.Mellapak250Y (                                                                                                                                                                                                        d = 8.65, H = 10), redeclare package MediumVapour =
            Media.H2O_CO2_Vap,                                                                                                                                                                                                        redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        n_elements = 10, redeclare model PressureLoss =
            PressureLoss.StructuredPackedColumn.NominalLinear (                                                                                                                                                                                                        Vdot_nom = 0.4, deltaP_nom(displayUnit = "Pa") = 3), redeclare model HeatTransferWall =
            Wall.ConstAlpha,                                                                                                                                                                                                        p_v_start_inlet = 100000, p_v_start_outlet = 100000, T_vap_start_bottom = 313.15, T_vap_start_top = 313.15, T_liq_start_bottom = 313.15, T_liq_start_top = 313.15, T_vapour_start = 323.15, T_liquid_start = 323.15, redeclare model ThermoEquilibrium =
            PhaseEquilibrium.H2O_CO2_MEA_startUp,                                                                                                                                                                                                        redeclare model FilmModel =
            FilmModel.StructuredPackedColumn.TrueEquilibriumStartUpCCSDesorption (                                                                                                                                   faktor_Ndot_v = 30000, redeclare model StateSelection =
                ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.StateSelection3))                                                                                                       annotation (
        Placement(transformation(extent = {{92, -55}, {182, 40}})));
      ThermalSeparation.Components.Pumps.idealPumpControlledVdot idealPumpControlledVdot(V_flow_start = 0.9, T = 1, dT = 0.5, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq)                                                                                                                                                    annotation (
        Placement(transformation(extent = {{19.5, -18.5}, {-19.5, 18.5}}, rotation = 90, origin = {61.5, -143.5})));
      ThermalSeparation.Components.LiquidVolumes.Tank tank(inertLiquid = {false, false, true}, d_volume = 8.65, level_start = 7.5, x_l_start = {0.8525, 0.0401, 0.1074}, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                               p_gas = 100000, T_start = 313.15) annotation (
        Placement(transformation(extent = {{8.5, -9}, {-8.5, 9}}, rotation = 0, origin = {136.5, -139})));
      ThermalSeparation.Components.LiquidVolumes.Tank tank1(inertLiquid = {false, false, true}, d_volume = 14.56, level_start = 4, x_l_start = {0.8525, 0.0401, 0.1074}, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                               p_gas = 100000, T_start = 308.15) annotation (
        Placement(transformation(extent = {{-109, -116}, {-76, -82}})));
      ThermalSeparation.Components.Pumps.idealPumpControlledVdot idealPumpControlledVdot1(V_flow_start = 0.9, T = 1, dT = 0.5, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq)                                                                                                                                                     annotation (
        Placement(transformation(extent = {{-23, -21}, {23, 21}}, rotation = 270, origin = {-31, -147})));
      Modelica.Blocks.Sources.RealExpression level(y = 4) annotation (
        Placement(transformation(extent = {{-90, -140}, {-68, -122}})));
      ThermalSeparation.Components.Condenser.FlashCondenser_CO2_H2O flash_Condenser_simple(T_out(displayUnit = "degC") = 313.15) annotation (
        Placement(transformation(extent = {{118, 82}, {170, 130}})));
      ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink3 annotation (
        Placement(transformation(extent = {{196, -12}, {208, 0}})));
      ThermalSeparation.Components.SourcesSinks.CombLiquid_x combLiquid_x(x_start = {0.8593, 0.0325, 0.1082}, redeclare package Medium =
            Media.H2O_CO2_MEA_Liq,                                                                                                                              T_start = 393.15) annotation (
        Placement(transformation(extent = {{-5, -5}, {5, 5}}, rotation = 270, origin = {-94, -67})));
      ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid_x(flowOption = Components.SourcesSinks.Enumerations.FlowOption.FlowMdot, Flow = 1021, use_Flow = true, x = {1, 0, 0}, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        T = 313.15) annotation (
        Placement(transformation(extent = {{-14, 118}, {-34, 98}})));
      Modelica.Blocks.Sources.RealExpression Water_in(y = -water_diff) annotation (
        Placement(transformation(extent = {{9.5, -9}, {-9.5, 9}}, rotation = 0, origin = {6.5, 103})));
      Modelica.Blocks.Sources.CombiTimeTable Q_Reboiler(smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint, offset = {1}, startTime = 0, table = [0, 0.0; 3000, 0; 3240, 2.13584e8; 5000, 2.13584e8; 20000, 2.13584e8; 23000, 1]) annotation (
        Placement(transformation(extent = {{-15, -15}, {15, 15}}, rotation = 180, origin = {217, -99})));
      ThermalSeparation.Components.Reboiler.KettleReboilerEq_StartUpCCS_dummy_p sumpV(m_Reb = 1000, eps_liq_init = 0.001, initOption_TepsXfixed = false, initOption_withoutT = false, initOption_standalone = false, initOption_startup_RebplusDes = true, x_init = 0.8593, inert_liq = {false, false, true}, nS = 2, startup = true, A = 50, H = 2, redeclare model InnerHT =
            HeatAndMassTransfer.HTResistance.NoHTResistance,                                                                                                                                                                                                        Vdot_nom = 100, beta_eps = 0.001, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare package MediumVapour =
            Media.H2O_CO2_Vap,                                                                                                                                                                                                        redeclare model Reaction =
            Reaction.ReactionEquilibrium.CO2_MEA,                                                                                                                                                                                                        init_option = ThermalSeparation.Components.Reboiler.InitOptionEq.initOption_startup_inert_Dyn, deltap_nom = 10, redeclare record Geometry =
            Geometry.StructuredPackedColumn.Mellapak250Y,                                                                                                                                                                                                        redeclare model ThermoEquilibrium =
            PhaseEquilibrium.CO2_CO2_MEA_StartUpReboiler_newFormulation,                                                                                                                                                                                                        T_init = 313.15, p_init = 100000, p_friggel = 25000) annotation (
        Placement(transformation(extent = {{126, -112}, {146, -92}})));
      Modelica.Blocks.Continuous.LimPID PIDLim(controllerType = Modelica.Blocks.Types.SimpleController.PI, yMax = 1.5, yMin = 0, k = 100) annotation (
        Placement(transformation(extent = {{-64, -102}, {-44, -82}})));
      ThermalSeparation.Components.SourcesSinks.CombLiquid_x combLiquid_x2(x_start = {0.8593, 0.0325, 0.1082}, redeclare package Medium =
            Media.H2O_CO2_MEA_Liq,                                                                                                                               T_start = 308.15) annotation (
        Placement(transformation(extent = {{-5, -5}, {5, 5}}, rotation = 180, origin = {-52, 87})));
      Modelica.Blocks.Math.IntegerToBoolean integerToBoolean annotation (
        Placement(transformation(extent = {{-224, -30}, {-204, -10}})));
      Modelica.Blocks.Math.IntegerToBoolean integerToBoolean1 annotation (
        Placement(transformation(extent = {{54, -46}, {74, -26}})));
      Modelica.Blocks.Sources.Ramp Vdot(duration = 240, offset = 0, height = 308.63, startTime = 3000) annotation (
        Placement(transformation(extent = {{-121, -76}, {-129, -68}})));
      Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(table = [0, 0.129, 0.026, 0.141, 0.704; 20000, 0.129, 0.026, 0.141, 0.704], extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint) annotation (
        Placement(transformation(extent = {{-198, -126}, {-178, -106}})));
      Modelica.Blocks.Sources.Ramp T_RG(height = -4.9, offset = 313.12, startTime = 20000, duration = 2550) annotation (
        Placement(transformation(extent = {{6, -6}, {-6, 6}}, rotation = 180, origin = {-190, -84})));
      ThermalSeparation.Components.SourcesSinks.SourceGas sourceGas_Vdot(flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowMdot, x = {0.0746, 0.0276, 0.1497, 0.7481}, use_T_in = true, use_x_in = true, redeclare package Medium =
            Media.H2O_O2_CO2_N2_Vap,                                                                                                                                                                                                        T = 312.15) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-163, -62})));
      Modelica.Blocks.Sources.IntegerExpression integerExpression annotation (
        Placement(transformation(extent = {{-258, -28}, {-238, -8}})));
      Modelica.Blocks.Sources.IntegerExpression integerExpression1 annotation (
        Placement(transformation(extent = {{4, -56}, {24, -36}})));
      Modelica.Blocks.Sources.BooleanExpression booleanExpression annotation (
        Placement(transformation(extent = {{-232, -2}, {-212, 18}})));
      Modelica.Blocks.Sources.BooleanExpression booleanExpression1 annotation (
        Placement(transformation(extent = {{50, -90}, {70, -70}})));
      inner SystemTS systemTS annotation (
        Placement(transformation(extent = {{-200, 160}, {-180, 180}})));
      ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink4 annotation (
        Placement(transformation(extent = {{184, 102}, {196, 114}})));
      Modelica.Blocks.Sources.CombiTimeTable Vdot1(smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint, startTime = 0, offset = {0}, table = [0, 0; 50, 0; 150, 0.929; 1500, 0.929; 20000, 0.929; 23000, 0; 26000, 0]) annotation (
        Placement(transformation(extent = {{-15, -15}, {15, 15}}, rotation = 180, origin = {101, -115})));
      Modelica.Blocks.Sources.CombiTimeTable mdot_in(smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint, startTime = 0, offset = {0}, table = [0, 0.0; 0, 0; 3000, 0; 3240, 308.63; 20000, 308.63; 23000, 308.63]) annotation (
        Placement(transformation(extent = {{-15, -15}, {15, 15}}, rotation = 180, origin = {-133, -123})));
      Modelica.Blocks.Sources.BooleanStep booleanStep(startValue = false, startTime = 21000) annotation (
        Placement(transformation(extent = {{54, -10}, {74, 10}})));
      Modelica.Blocks.Sources.BooleanStep booleanStep1(startTime = 23000, startValue = false) annotation (
        Placement(transformation(extent = {{-238, 24}, {-218, 44}})));
    equation
      connect(Water_in.y, sourceLiquid_x.Flow_In) annotation (
        Line(points = {{-3.95, 103}, {-14, 103}, {-14, 104}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(cooler.heatPort, ambientHeatSink1.heatPort1) annotation (
        Line(points = {{-77, 98.74}, {-77, 105.36}, {-76, 105.36}}, color = {188, 51, 69}, thickness = 1, smooth = Smooth.None));
      connect(Absorber.downStreamIn, cooler.coldLiquidOut) annotation (
        Line(points = {{-94.7, 65.1}, {-94.7, 86}, {-79.88, 86}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(tank1.portOut, idealPumpControlledVdot1.liquidIn) annotation (
        Line(points = {{-92.5, -115.32}, {-92.5, -147}, {-50.74, -147}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(idealPumpControlledVdot1.liquidOut, counterFlowHeatExchanger.coldLiquidIn) annotation (
        Line(points = {{-9.58, -147}, {-10, -147}, {-10, -44}, {-24, -44}, {-24, -12.18}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(counterFlowHeatExchanger.hotLiquidIn, idealPumpControlledVdot.liquidOut) annotation (
        Line(points = {{36, -12.72}, {36, -143.5}, {42.63, -143.5}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(tank.portOut, idealPumpControlledVdot.liquidIn) annotation (
        Line(points = {{136.5, -147.64}, {136.5, -166}, {106, -166}, {106, -143.5}, {78.89, -143.5}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(Desorber.upStreamOut, flash_Condenser_simple.gasPortIn) annotation (
        Line(points = {{105.5, 35.25}, {88, 35.25}, {88, 106}, {120.08, 106}}, color = {255, 127, 39}, thickness = 1, smooth = Smooth.None));
      connect(Absorber.heatPort, ambientHeatSink.heatPort1) annotation (
        Line(points = {{-87.84, 21}, {-72.96, 21}, {-72.96, 20}}, color = {188, 51, 69}, thickness = 1, smooth = Smooth.None));
      connect(Q_Reboiler.y[1], prescribedHeatFlow.Q_flow) annotation (
        Line(points = {{200.5, -99}, {182, -99}, {182, -102}, {174, -102}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(sumpV.gasPortOut, Desorber.upStreamIn) annotation (
        Line(points = {{130.2, -93}, {100, -93}, {100, -50.25}, {105.5, -50.25}}, color = {255, 127, 39}, thickness = 1, smooth = Smooth.None));
      connect(sumpV.liquidPortIn, Desorber.downStreamOut) annotation (
        Line(points = {{141, -93}, {166, -93}, {166, -50.25}, {168.5, -50.25}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(prescribedHeatFlow.port, sumpV.heatPort) annotation (
        Line(points = {{160, -102}, {144.2, -102}}, color = {191, 0, 0}, smooth = Smooth.None));
      connect(sinkGas.gasPortIn, Absorber.upStreamOut) annotation (
        Line(points = {{-162.5, 140.68}, {-162.5, 65.1}, {-163.3, 65.1}}, color = {255, 127, 39}, thickness = 1, smooth = Smooth.None));
      connect(counterFlowHeatExchanger.coldLiquidOut, Desorber.downStreamIn) annotation (
        Line(points = {{36, 5.64}, {36, 5.64}, {36, 66}, {170, 66}, {170, 42}, {170, 35.25}, {168.5, 35.25}}, color = {153, 217, 234}, thickness = 1));
      connect(combLiquid_x.liquidPortOut, tank1.portIn) annotation (
        Line(points = {{-94, -72}, {-94, -82.68}, {-92.5, -82.68}}, color = {153, 217, 234}, thickness = 1));
      connect(Absorber.downStreamOut, combLiquid_x.liquidPortIn2) annotation (
        Line(points = {{-94.7, -23.1}, {-94.7, -43.85}, {-96.5, -43.85}, {-96.5, -62}}, color = {153, 217, 234}, thickness = 1));
      connect(combLiquid_x2.liquidPortOut, cooler.hotLiquidIn) annotation (
        Line(points = {{-57, 87}, {-60.5, 87}, {-60.5, 86}, {-73.94, 86}}, color = {153, 217, 234}, thickness = 1));
      connect(combLiquid_x2.liquidPortIn1, counterFlowHeatExchanger.hotLiquidOut) annotation (
        Line(points = {{-47, 84.5}, {-47, 13.25}, {-24, 13.25}, {-24, 5.64}}, color = {153, 217, 234}, thickness = 1));
      connect(sourceLiquid_x.liquidPortOut, combLiquid_x2.liquidPortIn2) annotation (
        Line(points = {{-35.4, 108}, {-40, 108}, {-40, 89.5}, {-47, 89.5}}, color = {153, 217, 234}, thickness = 1));
      connect(sumpV.liquidPortOut, tank.portIn) annotation (
        Line(points = {{136, -111}, {136.5, -111}, {136.5, -130.36}}, color = {153, 217, 234}, thickness = 1));
      connect(tank1.y, PIDLim.u_s) annotation (
        Line(points = {{-77.65, -93.9}, {-78, -93.9}, {-78, -92}, {-66, -92}}, color = {0, 0, 127}));
      connect(PIDLim.y, idealPumpControlledVdot1.VdotInput) annotation (
        Line(points = {{-43, -92}, {-31, -92}, {-31, -119.4}}, color = {0, 0, 127}));
      connect(level.y, PIDLim.u_m) annotation (
        Line(points = {{-66.9, -131}, {-54, -131}, {-54, -104}}, color = {0, 0, 127}));
      connect(flash_Condenser_simple.liquidPortOut, combLiquid_x.liquidPortIn1) annotation (
        Line(points = {{142.7, 85.36}, {142.7, 72}, {-56, 72}, {-56, -62}, {-91.5, -62}}, color = {153, 217, 234}, thickness = 1));
      connect(ambientHeatSink3.heatPort1, Desorber.heatPort) annotation (
        Line(points = {{195.04, -6}, {174.8, -6}, {174.8, -7.5}}, color = {188, 51, 69}));
      connect(flash_Condenser_simple.gasPortOut, sinkGas1.gasPortIn) annotation (
        Line(points = {{142.96, 126.64}, {142.96, 138.32}, {146.5, 138.32}, {146.5, 148.64}}, color = {255, 127, 39}, thickness = 1));
      connect(integerToBoolean.y, Absorber.StartUp_signal) annotation (
        Line(points = {{-203, -20}, {-192, -20}, {-192, -8.89}, {-176.53, -8.89}}, color = {255, 0, 255}));
      connect(integerToBoolean1.y, Desorber.StartUp_signal) annotation (
        Line(points = {{75, -36}, {93.35, -36}, {93.35, -36.475}}, color = {255, 0, 255}));
      connect(sourceGas_Vdot.gasPortOut, Absorber.upStreamIn) annotation (
        Line(points = {{-163, -50.6}, {-163, -39.3}, {-163.3, -39.3}, {-163.3, -23.1}}, color = {255, 127, 39}, thickness = 1));
      connect(T_RG.y, sourceGas_Vdot.T_in) annotation (
        Line(points = {{-183.4, -84}, {-163, -84}, {-163, -72}}, color = {0, 0, 127}));
      connect(combiTimeTable.y, sourceGas_Vdot.x_in) annotation (
        Line(points = {{-177, -116}, {-159, -116}, {-159, -72}}, color = {0, 0, 127}));
      connect(integerExpression.y, integerToBoolean.u) annotation (
        Line(points = {{-237, -18}, {-232, -18}, {-232, -20}, {-226, -20}}, color = {255, 127, 0}));
      connect(integerExpression1.y, integerToBoolean1.u) annotation (
        Line(points = {{25, -46}, {25, -36}, {52, -36}}, color = {255, 127, 0}));
      connect(ambientHeatSink4.heatPort1, flash_Condenser_simple.heatPort) annotation (
        Line(points = {{183.04, 108}, {178, 108}, {178, 106}, {164.8, 106}}, color = {188, 51, 69}));
      connect(Vdot1.y[1], idealPumpControlledVdot.VdotInput) annotation (
        Line(points = {{84.5, -115}, {74.25, -115}, {74.25, -120.1}, {61.5, -120.1}}, color = {0, 0, 127}));
      connect(mdot_in.y[1], sourceGas_Vdot.Flow_in) annotation (
        Line(points = {{-149.5, -123}, {-167, -123}, {-167, -72}}, color = {0, 0, 127}));
      connect(booleanExpression.y, Absorber.ShutDown_signal) annotation (
        Line(points = {{-211, 8}, {-192, 8}, {-192, -1.05}, {-176.53, -1.05}}, color = {255, 0, 255}));
      connect(booleanExpression1.y, Desorber.ShutDown_signal) annotation (
        Line(points = {{71, -80}, {80, -80}, {80, -28.875}, {93.35, -28.875}}, color = {255, 0, 255}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -180}, {240, 180}})),
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -180}, {240, 180}})),
        experiment(StopTime = 20000, __Dymola_NumberOfIntervals = 5000));
    end Cycle_Heyden_StartUp;
  end Complex;

  package Pump
    model PumpSystem
      ThermalSeparation.Components.Pumps.idealPumpControlledVdot idealPumpControlledVdot(Vdot_l(start = 1), redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                  V_flow_start = 1, T_l_in(start = 308.15), T_l(start = 308.15), p(start = 130000)) annotation (
        Placement(transformation(extent = {{6, 6}, {-6, -6}}, rotation = 180, origin = {-6, 0})));
      ThermalSeparation.Components.SourcesSinks.SinkLiquid sinkLiquid(redeclare package Medium =
            Media.H2O_CO2_MEA_Liq)                                                                                      annotation (
        Placement(transformation(extent = {{20, 28}, {40, 48}})));
      ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid(Flow = 2, use_Flow = false, flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowVdot, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        T = 313.15, x = {0.875, 0.015, 0.11}) annotation (
        Placement(transformation(extent = {{-82, 6}, {-62, 26}})));
      Modelica.Blocks.Sources.Ramp ramp(height = 1, duration = 500, startTime = 500, offset = 1) annotation (
        Placement(transformation(extent = {{-42, 30}, {-32, 40}})));
      inner ThermalSeparation.SystemTS systemTS annotation (
        Placement(transformation(extent = {{-100, 80}, {-80, 100}})));
      ThermalSeparation.Components.LiquidVolumes.Tank tank(d_volume = 10, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                inertLiquid = {false, false, true}, x_l_start = {0.875, 0.015, 0.11}, level_start = 2, avoid_sum_mole = false, T_start = 313.15) annotation (
        Placement(transformation(extent = {{-40, -20}, {-20, 0}})));
    equation
      connect(sinkLiquid.liquidPortIn, idealPumpControlledVdot.liquidOut) annotation (
        Line(points = {{20.4, 38}, {12, 38}, {12, 6.12}, {-6, 6.12}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(ramp.y, idealPumpControlledVdot.VdotInput) annotation (
        Line(points = {{-31.5, 35}, {-31.5, 36}, {-13.2, 36}, {-13.2, 0}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(tank.portOut, idealPumpControlledVdot.liquidIn) annotation (
        Line(points = {{-30, -19.6}, {-30, -38}, {-6, -38}, {-6, -5.64}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(sourceLiquid.liquidPortOut, tank.portIn) annotation (
        Line(points = {{-60.6, 16}, {-46, 16}, {-46, 14}, {-30, 14}, {-30, -0.4}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}})),
        experiment(StopTime = 1500),
        __Dymola_experimentSetupOutput);
    end PumpSystem;
  end Pump;

  package Compressor
    model CompressorSimpleSystem
      Modelica.Blocks.Sources.Ramp ramp(duration = 500, startTime = 500, offset = 0.5, height = 0) annotation (
        Placement(transformation(extent = {{-62, -20}, {-52, -10}})));
      inner ThermalSeparation.SystemTS systemTS annotation (
        Placement(transformation(extent = {{-100, 80}, {-80, 100}})));
      ThermalSeparation.Components.Compressors.CompressorSimple compressorSimple(redeclare package Medium =
            ThermalSeparation.Media.IdealGasMixtures.N2_H2O,                                                                                                 useP = false, P_drive_const = 200000) annotation (
        Placement(transformation(extent = {{0, -4}, {20, 16}})));
      ThermalSeparation.Components.SourcesSinks.SinkGas sinkGas(redeclare package Medium =
            ThermalSeparation.Media.IdealGasMixtures.N2_H2O,                                                                                use_p = true, p = 500000) annotation (
        Placement(transformation(extent = {{34, 30}, {54, 50}})));
      ThermalSeparation.Components.SourcesSinks.SourceGas sourceGas(redeclare package Medium =
            ThermalSeparation.Media.IdealGasMixtures.N2_H2O,                                                                                    x = {0.95, 0.05}, T = 300, flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowVdot) annotation (
        Placement(transformation(extent = {{-36, -28}, {-16, -8}})));
    equation
      connect(sinkGas.gasPortIn, compressorSimple.gasPortOut) annotation (
        Line(points = {{34.4, 40}, {22, 40}, {22, 15.4}, {10, 15.4}}, color = {255, 127, 39}, thickness = 1, smooth = Smooth.None));
      connect(sourceGas.gasPortOut, compressorSimple.gasPortIn) annotation (
        Line(points = {{-14.6, -18}, {-2, -18}, {-2, -3.4}, {10, -3.4}}, color = {255, 127, 39}, thickness = 1, smooth = Smooth.None));
      connect(ramp.y, sourceGas.Flow_in) annotation (
        Line(points = {{-51.5, -15}, {-43.75, -15}, {-43.75, -14}, {-36, -14}}, color = {0, 0, 127}, smooth = Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics),
        experiment(StopTime = 1500),
        __Dymola_experimentSetupOutput);
    end CompressorSimpleSystem;
  end Compressor;

  package HeatExchanger
    model HeatExchangerSimple
      ThermalSeparation.Components.HeatExchanger.CcounterFlowHeatExchanger ccounterFlowHeatExchanger(c_l_start_hot = {5, 5, 5}, c_l_start_cold = {5, 5, 5}, ss_hot_c2 = false, ss_cold_c2 = false, surfaceArea = 3, alphaHot_nom_param = 3600, alphaCold_nom_param = 3600, redeclare package MediumLiquidCold =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare package MediumLiquidHot =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        ThotLiquidIn_start = 313.15, TcoldLiquidIn_start = 298.15, ThotLiquid_start = 302.15, TcoldLiquid_start = 308.15) annotation (
        Placement(transformation(extent = {{-12, -8}, {8, 12}})));
      ThermalSeparation.Components.SourcesSinks.SinkLiquid sinkLiquid(redeclare package Medium =
            Media.H2O_CO2_MEA_Liq)                                                                                      annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-42, -4})));
      ThermalSeparation.Components.SourcesSinks.SinkLiquid sinkLiquid1(redeclare package Medium =
            Media.H2O_CO2_MEA_Liq)                                                                                       annotation (
        Placement(transformation(extent = {{34, -14}, {54, 6}})));
      ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid1(Flow = 0.1, x = {0.87, 0.03, 0.1}, use_Flow = false, useT_In = true, flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowVdot, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        T = 363.15) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {32, 28})));
      ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid2(Flow = 0.1, x = {0.85, 0.05, 0.1}, useT_In = true, use_Flow = false, flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowVdot, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        T = 298.15) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-30, 28})));
      inner ThermalSeparation.SystemTS systemTS annotation (
        Placement(transformation(extent = {{-100, 80}, {-80, 100}})));
      Modelica.Blocks.Sources.Ramp ramp(duration = 500, startTime = 500, height = 15, offset = 25 + 273.15) annotation (
        Placement(transformation(extent = {{-86, 24}, {-66, 44}})));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 500, startTime = 500, height = 50, offset = 40 + 273.15) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {64, 30})));
    equation
      connect(sourceLiquid1.liquidPortOut, ccounterFlowHeatExchanger.hotLiquidIn) annotation (
        Line(points = {{20.6, 28}, {14, 28}, {14, 5.6}, {8, 5.6}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(sinkLiquid1.liquidPortIn, ccounterFlowHeatExchanger.coldLiquidOut) annotation (
        Line(points = {{34.4, -4}, {22, -4}, {22, -1.2}, {8, -1.2}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(sinkLiquid.liquidPortIn, ccounterFlowHeatExchanger.hotLiquidOut) annotation (
        Line(points = {{-32.4, -4}, {-22.2, -4}, {-22.2, -1.2}, {-12, -1.2}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(sourceLiquid2.liquidPortOut, ccounterFlowHeatExchanger.coldLiquidIn) annotation (
        Line(points = {{-18.6, 28}, {-16, 28}, {-16, 5.4}, {-12, 5.4}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(ramp.y, sourceLiquid2.T_In) annotation (
        Line(points = {{-65, 34}, {-52, 34}, {-52, 28}, {-40, 28}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(ramp1.y, sourceLiquid1.T_In) annotation (
        Line(points = {{53, 30}, {48, 30}, {48, 28}, {42, 28}}, color = {0, 0, 127}, smooth = Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics),
        experiment(StopTime = 1500),
        __Dymola_experimentSetupOutput);
    end HeatExchangerSimple;
  end HeatExchanger;

  package Reboiler

    model Reboiler_StartUp
      Components.SourcesSinks.SinkGas sinkGas(redeclare package Medium =
            Media.H2O_CO2_Vap,                                                              p = 200000) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-30, 30})));
      Components.SourcesSinks.SourceLiquid sourceLiquid(flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowVdot, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                       x = {0.87, 0.05, 0.1}, Flow = 1, T = 313.15, useT_In = true) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {32, 30})));
      Components.SourcesSinks.SinkLiquid sinkLiquid(redeclare package Medium =
            Media.H2O_CO2_MEA_Liq)                                                                    annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {0, -34})));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation (
        Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 180, origin = {27, 2})));
      inner SystemTS systemTS annotation (
        Placement(transformation(extent = {{-100, 80}, {-80, 100}})));
      Modelica.Blocks.Sources.CombiTimeTable Q_Reboiler(smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint, offset = {1}, startTime = 0, table = [0, 0.0; 3000, 0; 3240, 2.13584e8; 5000, 2.13584e8; 20000, 2.13584e8]) annotation (
        Placement(transformation(extent = {{-15, -15}, {15, 15}}, rotation = 180, origin = {69, 3})));
      Components.Reboiler.KettleReboilerEq_StartUpCCS_dummy_p sumpV1(m_Reb = 1000, initOption_TepsXfixed = false, initOption_withoutT = false, initOption_standalone = false, x_init = 0.8593, inert_liq = {false, false, true}, startup = true, A = 50, H = 2, redeclare model InnerHT =
            HeatAndMassTransfer.HTResistance.NoHTResistance,                                                                                                                                                                                                        Vdot_nom = 100, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare model Reaction =
            Reaction.ReactionEquilibrium.CO2_MEA,                                                                                                                                                                                                        init_option = ThermalSeparation.Components.Reboiler.InitOptionEq.initOption_startup_inert_Dyn, deltap_nom = 10, redeclare record Geometry =
            Geometry.StructuredPackedColumn.Mellapak250Y,                                                                                                                                                                                                        beta_eps = 0.001, initOption_startup_RebplusDes = true, initOption_startup_RebplusDesalone = false, mapping = [1, 1; 2, 2], nS = 2, redeclare package MediumVapour =
            Media.H2O_CO2_Vap,                                                                                                                                                                                                        redeclare model ThermoEquilibrium =
            PhaseEquilibrium.CO2_CO2_MEA_StartUpReboiler_newFormulation,                                                                                                                                                                                                        eps_liq_init = 0.0001, T_init = 313.15, p_init = 100000, p_friggel = 25000, p_dummy = 200000) annotation (
        Placement(transformation(extent = {{-8, -10}, {12, 10}})));
      Modelica.Blocks.Sources.Ramp T_solvent(height = 60, duration = 2000, offset = 313.15, startTime = 3240) annotation (
        Placement(transformation(extent = {{6, -6}, {-6, 6}}, rotation = 0, origin = {64, 36})));
    equation
      connect(prescribedHeatFlow.Q_flow, Q_Reboiler.y[1]) annotation (
        Line(points = {{34, 2}, {44, 2}, {44, 3}, {52.5, 3}, {52.5, 3}}, color = {0, 0, 127}));
      connect(sumpV1.liquidPortOut, sinkLiquid.liquidPortIn) annotation (
        Line(points = {{2, -9}, {2, -24.4}, {0, -24.4}}, color = {153, 217, 234}, thickness = 1));
      connect(sumpV1.gasPortOut, sinkGas.gasPortIn) annotation (
        Line(points = {{-3.8, 9}, {-3.8, 29.5}, {-20.4, 29.5}, {-20.4, 30}}, color = {255, 127, 39}, thickness = 1));
      connect(sourceLiquid.liquidPortOut, sumpV1.liquidPortIn) annotation (
        Line(points = {{20.6, 30}, {20.6, 29}, {7, 29}, {7, 9}}, color = {153, 217, 234}, thickness = 1));
      connect(prescribedHeatFlow.port, sumpV1.heatPort) annotation (
        Line(points = {{20, 2}, {16, 2}, {16, 0}, {10.2, 0}}, color = {191, 0, 0}));
      connect(sourceLiquid.T_In, T_solvent.y) annotation (
        Line(points = {{42, 30}, {50, 30}, {50, 36}, {57.4, 36}}, color = {0, 0, 127}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio = false)),
        Diagram(coordinateSystem(preserveAspectRatio = false)),
        experiment(StopTime = 25000, __Dymola_NumberOfIntervals = 5000));
    end Reboiler_StartUp;
  end Reboiler;

  package Column
    model ColumnSimple
      ThermalSeparation.Components.SourcesSinks.SinkGas sinkGas(redeclare package Medium =
            Media.H2O_CO2_Vap,                                                                                p = 200000) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-50, 72})));
      ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid(flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowVdot, Flow = 0.001472, x = {0.87, 0.03, 0.1}, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        T = 378.15) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {28, 62})));
      ThermalSeparation.Components.SourcesSinks.SinkLiquid sinkLiquid(redeclare package Medium =
            Media.H2O_CO2_MEA_Liq)                                                                                      annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {4, -34})));
      inner ThermalSeparation.SystemTS systemTS annotation (
        Placement(transformation(extent = {{-100, 80}, {-80, 100}})));
      ThermalSeparation.Components.Columns.StructuredPackedColumn Desorber( redeclare model HeatTransferWall =
            ThermalSeparation.Wall.Adiabatic,                                                                                                    redeclare model Holdup =
            ThermalSeparation.Holdup.StructuredPackedColumn.StichlmairStat,                                                                                                                                                               redeclare model PressureLoss =
            ThermalSeparation.PressureLoss.StructuredPackedColumn.NominalLinear,                                                                                                                                                                                                        redeclare package MediumVapour =
            Media.H2O_CO2_Vap,                                                                                                                                                                                                        redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare model ThermoEquilibrium =
            PhaseEquilibrium.H2O_CO2_MEA,                                                                                                                                                                                                        redeclare model Reaction =
            Reaction.ReactionEquilibrium.CO2_MEA,                                                                                                                                                                                                        redeclare record Geometry =
            Geometry.StructuredPackedColumn.Mellapak250Y (                                                                                                                                                                                                        d = 1, H = 1), redeclare model InitOption =
            Components.Columns.BaseClasses.Initialization.Init_EquilibriumFilm,                                                                                                                                                                                                        redeclare model BalanceEquations =
            BalanceEquations.StructuredPackedColumn.NonEquilibrium.TwoPhaseVarState (                                                                                                                                redeclare model FilmModel =
                ThermalSeparation.FilmModel.StructuredPackedColumn.TrueEquilibrium),                                                                                                                                                                                                        T_l_profile = false, T_liquid_start = 387.15, T_v_profile = false, T_vap_start_bottom = 383.15, T_vap_start_top = 373.15, T_vapour_start = 387.15, eps_liq_start = 0.05, inertLiquid = {false, false, true},mapping = {{1, 1}, {2, 2}}, nS = 2, n_elements = 1, p_v_in(start = 200500), p_v_start_inlet = 200500, p_v_start_outlet = 200000, useHomotopy = false, wettedInitial = true, x_l_profile = false, x_l_start_const = {0.87, 0.03, 0.1}, x_l_start_in = {0.87, 0.03, 0.1}, x_l_start_out = {0.87, 0.03, 0.1}, x_v_start_const = {0.6, 0.4}) annotation (
        Placement(transformation(extent = {{-36, -9}, {10, 37}})));
      ThermalSeparation.Components.SourcesSinks.SourceGas sourceGas(flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowNdot, Flow = 7, x = {0.7, 0.3}, use_Flow = false, redeclare package Medium =
            Media.H2O_CO2_Vap,                                                                                                                                                                                                        T = 393.15) annotation (
        Placement(transformation(extent = {{-60, -38}, {-40, -18}})));
      ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink annotation (
        Placement(transformation(extent = {{20, 2}, {32, 14}})));
    equation
      connect(sourceLiquid.liquidPortOut, Desorber.downStreamIn) annotation (
        Line(points = {{16.6, 62}, {10, 62}, {10, 34.7}, {3.1, 34.7}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(sinkGas.gasPortIn, Desorber.upStreamOut) annotation (
        Line(points = {{-40.4, 72}, {-36, 72}, {-36, 34.7}, {-29.1, 34.7}}, color = {255, 127, 39}, thickness = 1, smooth = Smooth.None));
      connect(sourceGas.gasPortOut, Desorber.upStreamIn) annotation (
        Line(points = {{-38.6, -28}, {-34, -28}, {-34, -6.7}, {-29.1, -6.7}}, color = {255, 127, 39}, thickness = 1, smooth = Smooth.None));
      connect(Desorber.downStreamOut, sinkLiquid.liquidPortIn) annotation (
        Line(points = {{3.1, -6.7}, {3.1, -15.35}, {4, -15.35}, {4, -24.4}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(ambientHeatSink.heatPort1, Desorber.heatPort) annotation (
        Line(points = {{19.04, 8}, {14, 8}, {14, 14}, {6.32, 14}}, color = {188, 51, 69}, smooth = Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics),
        experiment(StopTime = 1500),
        __Dymola_experimentSetupOutput);
    end ColumnSimple;

    model AbsorptionColumnStartUp
      Components.Columns.StructuredPackedColumn_newStartUpShutDown Absorber(inertVapour = {false, true, false, true}, inertLiquid = {false, false, true}, nS = 2, mapping = {{1, 1}, {3, 2}}, wettedInitial = false, redeclare model Holdup =
            ThermalSeparation.Holdup.StructuredPackedColumn.StichlmairStat,                                                                                                                                                                                                        considerStartUp = true, StartUp_CCS = true, friggelfaktor = 0.2e5, smooth_startUp = false, delay_startUp = 750, lowBoilingPoint = {false, true, true, true}, redeclare model InitOption =
            ThermalSeparation.Components.Columns.BaseClasses.Initialization.DyncapStartUpAbsorption,                                                                                                                                                                                                        eps_liq_start = 0.001, x_v_start_const = {0.048, 0.2, 0.001, 0.751}, switchingCondition_Boiling = false, x_v_dummy = {1, 0, 0, 0}, y_PID = 1000, gain = 100, Vdot_startUp_pressure = 5, x_l_start_const = {0.8593, 0.0325, 0.1082}, redeclare record Geometry =
            ThermalSeparation.Geometry.StructuredPackedColumn.Mellapak250Y (                                                                                                                                                                                                        d = 14.56, H = 15), redeclare model PressureLoss =
            ThermalSeparation.PressureLoss.StructuredPackedColumn.NominalLinear (                                                                                                                                    deltaP_nom(displayUnit = "Pa") = 3), redeclare package MediumVapour =
            ThermalSeparation.Media.H2O_O2_CO2_N2_Vap,                                                                                                                                                                                                        redeclare package MediumLiquid =
            ThermalSeparation.Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare model ThermoEquilibrium =
            PhaseEquilibrium.H2O_CO2_MEA_startUp (                                                                                                                                                                                                        factor_K = {1, 1.03, 1}), redeclare model Reaction =
            Reaction.ReactionEquilibrium.CO2_MEA,                                                                                                                                                                                                        redeclare model FilmModel =
            FilmModel.StructuredPackedColumn.TrueEquilibriumStartUpCCSAbsorption (                                                                                                                                   faktor_Ndot_v = 30000), redeclare model HeatTransferWall =
            Wall.Adiabatic,
        n_elements=5,
        p_v_start_inlet=100000,
        p_v_start_outlet=100000,
        T_vapour_start=313.15,
        T_liquid_start=313.15,
        p_v_in(start=100000))                                                                                                                                                                                                         annotation (
        Placement(transformation(extent = {{-58, -44}, {40, 54}})));
      Modelica.Blocks.Math.IntegerToBoolean integerToBoolean annotation (
        Placement(transformation(extent = {{-102, -32}, {-82, -12}})));
      Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(table = [0, 0.129, 0.026, 0.141, 0.704; 20000, 0.129, 0.026, 0.141, 0.704], extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint) annotation (
        Placement(transformation(extent = {{-76, -128}, {-56, -108}})));
      Modelica.Blocks.Sources.Ramp T_RG(height = -4.9, offset = 313.12, startTime = 20000, duration = 2550) annotation (
        Placement(transformation(extent = {{6, -6}, {-6, 6}}, rotation = 180, origin = {-68, -86})));
      Modelica.Blocks.Sources.IntegerExpression integerExpression annotation (
        Placement(transformation(extent = {{-136, -30}, {-116, -10}})));
      Modelica.Blocks.Sources.BooleanExpression booleanExpression annotation (
        Placement(transformation(extent = {{-110, -4}, {-90, 16}})));
      Modelica.Blocks.Sources.BooleanStep booleanStep1(startTime = 23000, startValue = false) annotation (
        Placement(transformation(extent = {{-116, 22}, {-96, 42}})));
      Components.SourcesSinks.SourceGas sourceGas_Vdot(flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowMdot, x = {0.0746, 0.0276, 0.1497, 0.7481}, use_T_in = true, use_x_in = true, redeclare package Medium =
            Media.H2O_O2_CO2_N2_Vap,                                                                                                                                                                                                        T = 312.15) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-41, -64})));
      Modelica.Blocks.Sources.CombiTimeTable mdot_in(smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint, startTime = 0,
        table=[0,0.0; 0,0; 3000,0; 3240,308.63; 20000,308.63; 20100,308.63],
        offset={0.1})                                                                                                                                                                                                         annotation (
        Placement(transformation(extent = {{-15, -15}, {15, 15}}, rotation = 180, origin = {7, -85})));
      Components.SourcesSinks.SinkGas sinkGas(redeclare package Medium =
            Media.H2O_O2_CO2_N2_Vap,                                                              use_p = true,
        p=100000)                                                                                                           annotation (
        Placement(transformation(extent = {{-17, -21.5}, {17, 21.5}}, rotation = 90, origin = {-44.5, 93})));
      Components.SourcesSinks.SourceLiquid sourceLiquid_x(flowOption = Components.SourcesSinks.Enumerations.FlowOption.FlowMdot, Flow = 1021, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                    x = {0.85, 0.03, 0.12}, use_Flow = false, T = 313.15) annotation (
        Placement(transformation(extent = {{88, 74}, {68, 54}})));
      Components.SourcesSinks.SinkLiquid sinkLiquid(redeclare package Medium =
            Media.H2O_CO2_MEA_Liq)                                                                    annotation (
        Placement(transformation(extent = {{72, -64}, {92, -44}})));
      Components.SourcesSinks.AmbientHeatSink ambientHeatSink annotation (
        Placement(transformation(extent = {{66, 0}, {78, 12}})));
      inner SystemTS systemTS annotation (
        Placement(transformation(extent = {{-160, 120}, {-140, 140}})));
    equation
      connect(integerToBoolean.y, Absorber.StartUp_signal) annotation (
        Line(points = {{-81, -22}, {-70, -22}, {-70, -24.89}, {-56.53, -24.89}}, color = {255, 0, 255}));
      connect(T_RG.y, sourceGas_Vdot.T_in) annotation (
        Line(points = {{-61.4, -86}, {-41, -86}, {-41, -74}}, color = {0, 0, 127}));
      connect(integerExpression.y, integerToBoolean.u) annotation (
        Line(points = {{-115, -20}, {-110, -20}, {-110, -22}, {-104, -22}}, color = {255, 127, 0}));
      connect(booleanStep1.y, Absorber.ShutDown_signal) annotation (
        Line(points = {{-95, 32}, {-76, 32}, {-76, -17.05}, {-56.53, -17.05}}, color = {255, 0, 255}));
      connect(sourceGas_Vdot.gasPortOut, Absorber.upStreamIn) annotation (
        Line(points = {{-41, -52.6}, {-41, -44.3}, {-43.3, -44.3}, {-43.3, -39.1}}, color = {255, 127, 39}, thickness = 1));
      connect(sinkGas.gasPortIn, Absorber.upStreamOut) annotation (
        Line(points = {{-44.5, 76.68}, {-44.5, 62.34}, {-43.3, 62.34}, {-43.3, 49.1}}, color = {255, 127, 39}, thickness = 1));
      connect(mdot_in.y[1], sourceGas_Vdot.Flow_in) annotation (
        Line(points = {{-9.5, -85}, {-45, -85}, {-45, -74}}, color = {0, 0, 127}));
      connect(combiTimeTable.y, sourceGas_Vdot.x_in) annotation (
        Line(points = {{-55, -118}, {-46, -118}, {-37, -118}, {-37, -74}}, color = {0, 0, 127}));
      connect(sourceLiquid_x.liquidPortOut, Absorber.downStreamIn) annotation (
        Line(points = {{66.6, 64}, {25.3, 64}, {25.3, 49.1}}, color = {153, 217, 234}, thickness = 1));
      connect(sinkLiquid.liquidPortIn, Absorber.downStreamOut) annotation (
        Line(points = {{72.4, -54}, {50, -54}, {50, -52}, {25.3, -52}, {25.3, -39.1}}, color = {153, 217, 234}, thickness = 1));
      connect(ambientHeatSink.heatPort1, Absorber.heatPort) annotation (
        Line(points = {{65.04, 6}, {32.16, 6}, {32.16, 5}}, color = {188, 51, 69}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-160, -140}, {160, 140}})),
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-160, -140}, {160, 140}})),
        experiment(StopTime = 10000));
    end AbsorptionColumnStartUp;

    model DesorptionColumnStartUp
      Components.SourcesSinks.AmbientHeatSink ambientHeatSink2 annotation (
        Placement(transformation(extent = {{50, 3}, {62, 12}})));
      Components.SourcesSinks.SinkGas sinkGas1(redeclare package Medium =
            Media.H2O_CO2_Vap,                                                               p = 200000) annotation (
        Placement(transformation(extent = {{-16, -17.5}, {16, 17.5}}, rotation = 90, origin = {66.5, 152})));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation (
        Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 180, origin = {89, -122})));
      Components.Columns.StructuredPackedColumn_newStartUpShutDown Desorber(mapping = {{1, 1}, {2, 2}}, inertLiquid = {false, false, true}, T_l_profile = false, x_l_start_in = {0.8, 0.1, 0.1}, x_l_start_out = {0.5, 0.2, 0.3}, x_l_profile = false, T_v_profile = false, nS = 2, inertVapour = {false, false}, wettedInitial = false, eps_liq_start = 0.001, considerStartUp = true, StartUp_CCS = true, switchingCondition_Boiling = true, lowBoilingPoint = {false, true}, redeclare model InitOption =
            ThermalSeparation.Components.Columns.BaseClasses.Initialization.DyncapStartUpDesorption,                                                                                                                                                                                                        redeclare model Holdup =
            ThermalSeparation.Holdup.StructuredPackedColumn.StichlmairStat,                                                                                                                                                                                                        x_v_start_const = {0, 1}, y_PID = 1000, Vdot_startUp_pressure = 5, gain = 100, smooth_startUp = true, k = 0.04, delay_startUp = 100, friggelfaktor = 0.7e5, redeclare model Reaction =
            Reaction.ReactionEquilibrium.CO2_MEA,                                                                                                                                                                                                        redeclare record Geometry =
            Geometry.StructuredPackedColumn.Mellapak250Y (                                                                                                                                                                                                        d = 8.65, H = 10), redeclare package MediumVapour =
            Media.H2O_CO2_Vap,                                                                                                                                                                                                        redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        n_elements = 10, redeclare model PressureLoss =
            PressureLoss.StructuredPackedColumn.NominalLinear (                                                                                                                                                                                                        Vdot_nom = 0.4, deltaP_nom(displayUnit = "Pa") = 3), redeclare model HeatTransferWall =
            Wall.ConstAlpha,                                                                                                                                                                                                        redeclare model FilmModel =
            FilmModel.StructuredPackedColumn.TrueEquilibriumStartUpCCSDesorption (                                                                                                                                   faktor_Ndot_v = 30000, redeclare model StateSelection =
                ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.None),                                                                                                                                                                                                        redeclare model ThermoEquilibrium =
            PhaseEquilibrium.H2O_CO2_MEA_startUp,                                                                                                                                                                                                        p_v_start_inlet = 100000, p_v_start_outlet = 100000, x_l_start_const = {0.8593, 0.0325, 0.1082}, T_vap_start_bottom = 313.15, T_vap_start_top = 313.15, T_liq_start_bottom = 313.15, T_liq_start_top = 313.15, T_vapour_start = 323.15, T_liquid_start = 323.15) annotation (
        Placement(transformation(extent = {{12, -67}, {102, 28}})));
      Components.Condenser.FlashCondenser_CO2_H2O flash_Condenser_simple(T_out(displayUnit = "degC") = 313.15) annotation (
        Placement(transformation(extent = {{38, 70}, {90, 118}})));
      Components.SourcesSinks.AmbientHeatSink ambientHeatSink3 annotation (
        Placement(transformation(extent = {{116, -24}, {128, -12}})));
      Modelica.Blocks.Sources.CombiTimeTable Q_Reboiler(smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint, offset = {1}, startTime = 0, table = [0, 0.0; 3000, 0; 3240, 2.13584e8; 5000, 2.13584e8; 20000, 2.13584e8; 23000, 2.13584e8]) annotation (
        Placement(transformation(extent = {{-15, -15}, {15, 15}}, rotation = 180, origin = {137, -111})));
      Components.Reboiler.KettleReboilerEq_StartUpCCS_dummy_p sumpV(m_Reb = 1000, eps_liq_init = 0.001, initOption_TepsXfixed = false, initOption_withoutT = false, initOption_standalone = false, initOption_startup_RebplusDes = true, x_init = 0.8593, inert_liq = {false, false, true}, nS = 2, startup = true, A = 50, H = 2, redeclare model InnerHT =
            HeatAndMassTransfer.HTResistance.NoHTResistance,                                                                                                                                                                                                        Vdot_nom = 100, redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                                                                                                                                        redeclare package MediumVapour =
            Media.H2O_CO2_Vap,                                                                                                                                                                                                        redeclare model Reaction =
            Reaction.ReactionEquilibrium.CO2_MEA,                                                                                                                                                                                                        init_option = ThermalSeparation.Components.Reboiler.InitOptionEq.initOption_startup_inert_Dyn, deltap_nom = 10, redeclare record Geometry =
            Geometry.StructuredPackedColumn.Mellapak250Y,                                                                                                                                                                                                        redeclare model ThermoEquilibrium =
            PhaseEquilibrium.CO2_CO2_MEA_StartUpReboiler_newFormulation,                                                                                                                                                                                                        beta_eps = 0.0001, T_init = 313.15, p_init = 100000, p_friggel = 25000) annotation (
        Placement(transformation(extent = {{48, -132}, {68, -112}})));
      Components.SourcesSinks.AmbientHeatSink ambientHeatSink4 annotation (
        Placement(transformation(extent = {{104, 90}, {116, 102}})));
      Components.SourcesSinks.SinkLiquid sinkLiquid(redeclare package Medium =
            Media.H2O_CO2_MEA_Liq)                                                                    annotation (
        Placement(transformation(extent = {{64, -160}, {84, -140}})));
      Modelica.Blocks.Sources.BooleanExpression booleanExpression1 annotation (
        Placement(transformation(extent = {{-38, -34}, {-18, -14}})));
      Modelica.Blocks.Sources.IntegerExpression integerExpression1 annotation (
        Placement(transformation(extent = {{-84, -64}, {-64, -44}})));
      Modelica.Blocks.Math.IntegerToBoolean integerToBoolean1 annotation (
        Placement(transformation(extent = {{-36, -64}, {-16, -44}})));
      inner SystemTS systemTS annotation (
        Placement(transformation(extent = {{-100, 160}, {-80, 180}})));
      Components.SourcesSinks.SourceLiquid sourceLiquid_x(redeclare package MediumLiquid =
            Media.H2O_CO2_MEA_Liq,                                                                                flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowVdot, Flow = 1, use_Flow = false, useT_In = true, T = 313.15, x = {0.85, 0.03, 0.12}) annotation (
        Placement(transformation(extent = {{146, 74}, {126, 54}})));
      Components.SourcesSinks.CombLiquid_x combLiquid_x2(x_start = {0.8593, 0.0325, 0.1082}, redeclare package Medium =
            Media.H2O_CO2_MEA_Liq,                                                                                                             T_start = 308.15) annotation (
        Placement(transformation(extent = {{-5, -5}, {5, 5}}, rotation = 270, origin = {110, 43})));
      Modelica.Blocks.Sources.Ramp T_solvent(height = 60, duration = 2000, offset = 313.15, startTime = 3240) annotation (
        Placement(transformation(extent = {{6, -6}, {-6, 6}}, rotation = 0, origin = {182, 70})));
    equation
      connect(Desorber.upStreamOut, flash_Condenser_simple.gasPortIn) annotation (
        Line(points = {{25.5, 23.25}, {8, 23.25}, {8, 94}, {40.08, 94}}, color = {255, 127, 39}, thickness = 1, smooth = Smooth.None));
      connect(Q_Reboiler.y[1], prescribedHeatFlow.Q_flow) annotation (
        Line(points = {{120.5, -111}, {102, -111}, {102, -122}, {96, -122}}, color = {0, 0, 127}, smooth = Smooth.None));
      connect(sumpV.liquidPortIn, Desorber.downStreamOut) annotation (
        Line(points = {{63, -113}, {86, -113}, {86, -62.25}, {88.5, -62.25}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(prescribedHeatFlow.port, sumpV.heatPort) annotation (
        Line(points = {{82, -122}, {66.2, -122}}, color = {191, 0, 0}, smooth = Smooth.None));
      connect(ambientHeatSink3.heatPort1, Desorber.heatPort) annotation (
        Line(points = {{115.04, -18}, {94.8, -18}, {94.8, -19.5}}, color = {188, 51, 69}));
      connect(flash_Condenser_simple.gasPortOut, sinkGas1.gasPortIn) annotation (
        Line(points = {{62.96, 114.64}, {62.96, 126.32}, {66.5, 126.32}, {66.5, 136.64}}, color = {255, 127, 39}, thickness = 1));
      connect(ambientHeatSink4.heatPort1, flash_Condenser_simple.heatPort) annotation (
        Line(points = {{103.04, 96}, {98, 96}, {98, 94}, {84.8, 94}}, color = {188, 51, 69}));
      connect(sinkLiquid.liquidPortIn, sumpV.liquidPortOut) annotation (
        Line(points = {{64.4, -150}, {58, -150}, {58, -131}}, color = {153, 217, 234}, thickness = 1));
      connect(booleanExpression1.y, Desorber.ShutDown_signal) annotation (
        Line(points = {{-17, -24}, {-8, -24}, {-8, -40.875}, {13.35, -40.875}}, color = {255, 0, 255}));
      connect(integerExpression1.y, integerToBoolean1.u) annotation (
        Line(points = {{-63, -54}, {-38, -54}}, color = {255, 127, 0}));
      connect(integerToBoolean1.y, Desorber.StartUp_signal) annotation (
        Line(points = {{-15, -54}, {0, -54}, {0, -48.475}, {13.35, -48.475}}, color = {255, 0, 255}));
      connect(combLiquid_x2.liquidPortOut, Desorber.downStreamIn) annotation (
        Line(points = {{110, 38}, {88.5, 38}, {88.5, 23.25}}, color = {153, 217, 234}, thickness = 1));
      connect(combLiquid_x2.liquidPortIn2, flash_Condenser_simple.liquidPortOut) annotation (
        Line(points = {{107.5, 48}, {84.75, 48}, {84.75, 73.36}, {62.7, 73.36}}, color = {153, 217, 234}, thickness = 1));
      connect(combLiquid_x2.liquidPortIn1, sourceLiquid_x.liquidPortOut) annotation (
        Line(points = {{112.5, 48}, {120, 48}, {120, 64}, {124.6, 64}}, color = {153, 217, 234}, thickness = 1));
      connect(T_solvent.y, sourceLiquid_x.T_In) annotation (
        Line(points = {{175.4, 70}, {166, 70}, {166, 64}, {146, 64}}, color = {0, 0, 127}));
      connect(sumpV.gasPortOut, Desorber.upStreamIn) annotation (
        Line(points = {{52.2, -113}, {52.2, -88.5}, {25.5, -88.5}, {25.5, -62.25}}, color = {255, 127, 39}, thickness = 1));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -160}, {260, 180}})),
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -160}, {260, 180}})),
        experiment(StopTime = 25000));
    end DesorptionColumnStartUp;
  end Column;

  package LiquidCooler
    model LiquidCoolerSimple
      ThermalSeparation.Components.SourcesSinks.SinkLiquid sinkLiquid1(redeclare package Medium =
            ThermalSeparation.Media.WaterBasedLiquid.N2_H2O)                                                                                       annotation (
        Placement(transformation(extent = {{46, 14}, {66, 34}})));
      ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid2(flowOption = ThermalSeparation.Components.SourcesSinks.Enumerations.FlowOption.FlowVdot, Flow = 0.1, use_Flow = false, redeclare package MediumLiquid =
            ThermalSeparation.Media.WaterBasedLiquid.N2_H2O,                                                                                                                                                                                                        x = {0.1, 0.9}, useT_In = false, T = 353.15) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-30, 28})));
      inner ThermalSeparation.SystemTS systemTS annotation (
        Placement(transformation(extent = {{-100, 80}, {-80, 100}})));
      ThermalSeparation.Components.HeatExchanger.LiquidCooler liquidCooler(redeclare package MediumLiquid =
            ThermalSeparation.Media.WaterBasedLiquid.N2_H2O,                                                                                                 T_set = 298.15) annotation (
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {16, 38})));
      ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink annotation (
        Placement(transformation(extent = {{24, 56}, {38, 70}})));
    equation
      connect(sourceLiquid2.liquidPortOut, liquidCooler.hotLiquidIn) annotation (
        Line(points = {{-18.6, 28}, {12.6, 28}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(liquidCooler.coldLiquidOut, sinkLiquid1.liquidPortIn) annotation (
        Line(points = {{19.2, 28}, {32, 28}, {32, 24}, {46.4, 24}}, color = {153, 217, 234}, thickness = 1, smooth = Smooth.None));
      connect(ambientHeatSink.heatPort1, liquidCooler.heatPort) annotation (
        Line(points = {{22.88, 63}, {22.88, 52.5}, {16, 52.5}, {16, 46.2}}, color = {188, 51, 69}, smooth = Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics),
        experiment(StopTime = 1500),
        __Dymola_experimentSetupOutput);
    end LiquidCoolerSimple;
  end LiquidCooler;
end Testing;
