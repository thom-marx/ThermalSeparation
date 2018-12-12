within ThermalSeparation.Components.GasLiquidVolumes;
package Tests
  extends Icons.Library.Red;

  model TestLSF
  // Lauftest: 25.7.2011
  // Lauftest: 14.7.2011
    inner SystemTS systemTS 
      annotation (Placement(transformation(extent={{60,60},{80,80}})));
    LSF_deltaP leanSolventFlash(
      x_v_start={0.5,0.5},
      c_v_guess={50,50},
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
      initEQ=true,
      redeclare package MediumLiquid = 
          ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS,
      p_start=200000) 
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
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
      prescribedHeatFlow2 
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

  model TestLSF_noDeltaP
    // Lauftest: 25.7.2011
    //Lauftest: 14.7.2011
  //Lauftest: 1.7.2011
    inner SystemTS systemTS 
      annotation (Placement(transformation(extent={{60,60},{80,80}})));
    LSF_noDeltaP leanSolventFlash(
      x_v_start={0.5,0.5},
      c_v_guess={50,50},
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
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
      prescribedHeatFlow2 
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

  model Control3_2_tank_torricelli
    "Control Strategie 3 (idealPump1 regelt Abscheiderate, Qdot_reboiler wird gestellt, idealPump regelt Füllstand)"
    //Lauftest: 28.7.2011, 18000 s in 498 s (Kevin): geht erst, seit man an dem PI-Regler für den StartUp rumgespielt hat
    //Lauftest: 18.7.2011, 18000 s in 498 s
    //Lauftest: 26.5.2011, 18000 s in 700 s (Torricelli), wenn Umrechnung von x und c nicht mit v sondern mit MM/rho gemacht wird
    //Lauftest: 18.5.2011 / 20.5.2011 18000s in 600 s (Torricelli)
    //Lauftest: 11.5.2011, die ersten 6000 s werden in ca. 360 s gerechnet (Torricelli)
  //Lauftest: 21.4.2011
    parameter Integer switch = 8000
      "time in s where change is applied to boundary conditions";
          Modelica.SIunits.MassFlowRate mdot_desorber_water_out=Desorber.Vdot_v[
        Desorber.n]*Desorber.rho_v[Desorber.n]*Desorber.x_v[Desorber.n, 1]*0.018/
        Desorber.MM_v[Desorber.n];
      Modelica.SIunits.MassFlowRate mdot_desorber_CO2_out=Desorber.Vdot_v[
        Desorber.n]*Desorber.rho_v[Desorber.n]*Desorber.x_v[Desorber.n, 2]*0.044/
        Desorber.MM_v[Desorber.n];
      Modelica.SIunits.MassFlowRate mdot_water_out=mdot_desorber_CO2_out*(1/
        X_CO2_out - 1);
      Real X_CO2_out= 1-0.013;
      Modelica.SIunits.MassFlowRate mdot_desorber_water_in=
        mdot_desorber_water_out - mdot_water_out;

  ThermalSeparation.Components.SourcesSinks.SinkGas sinkGas1(redeclare package
        Medium = 
          ApplicationsThermalSeparation.Media.IdealGasMixtures.H2O_CO2_forPI,
                                                                         p=100000) 
                  annotation (Placement(transformation(extent={{112,158},{132,178}},
            rotation=0)));
    ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink2 
      annotation (Placement(transformation(extent={{188,116},{208,136}},
                                                                    rotation=0)));

    ThermalSeparation.Components.Pumps.idealPumpControlledVdot idealPump1(
      redeclare package MediumLiquid = 
          ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS_WaterProp,
      T_l(start=324.15),
      p(start=600000))                   annotation (Placement(transformation(
          extent={{12,-12},{-12,12}},
          rotation=90,
          origin={132,-54})));

    ThermalSeparation.Components.LiquidVolumes.Sump sump(
      d_volume=0.07,
      T_start(displayUnit="K") = 374.81,
      zeta=30,
      x_start={0.0397,0.9446,0.0157},
      redeclare package MediumLiquid = 
          ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS_WaterProp,
      p(start=105000),
      level_start=5) 
      annotation (Placement(transformation(extent={{166,70},{184,90}})));

    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
      prescribedHeatFlow1 
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={252,38})));
    Modelica.Blocks.Sources.Ramp ramp3(
      startTime=switch,
      duration=10,
      height=-0.3*25.3e3,
      offset=25.3e3) 
                   annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={288,38})));
    ThermalSeparation.Components.GasLiquidVolumes.Reboiler reboiler(
      x_v_start={0.1,0.9},
      c_v_guess={60,1},
      mapping={{1,2},{2,1}},
      inertVapour={false,false},
      zeta=8.6,
      eps_liq_start=0.15,
      x_l_start={0.0397,0.9446,0.0157},
      redeclare package MediumVapour = 
          ApplicationsThermalSeparation.Media.IdealGasMixtures.H2O_CO2_forPI,
      redeclare package MediumLiquid = 
          ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS_WaterProp,
      T_vapour_start(displayUnit="K") = 373,
      T_liquid_start(displayUnit="K") = 373,
      length_HX=0.07,
      redeclare model Reaction = 
          ApplicationsThermalSeparation.Obsolete.CO2_Siemens (
          moleBalance=true,
          a_R={0},
          b_R={-50e3}),
      p_out(start=101500),
      p_start=105500) 
      annotation (Placement(transformation(extent={{200,24},{230,58}})));

  ThermalSeparation.Components.Columns.StructuredPackedColumn Desorber(
      hasVapourFeed=false,
      nS=2,
      stageLiquidFeed={20},
      hasLiquidFeed=false,
      c_l_guess={500,50000,500},
      inertLiquid={false,false,true},
      c_v_guess={60,1},
      mapping={{1,2},{2,1}},
      inertVapour={false,false},
      wettedInitial=false,
      x_v_start_const={0.9,0.1},
      redeclare model Holdup = 
          ThermalSeparation.Holdup.StructuredPackedColumn.Rocha (
            hu_stat_const=false, frickel_dyn=1),
      redeclare model HeatTransferWall = ThermalSeparation.Wall.ConstQdot (stat=
              true, Qdot_set=-1000),
      x_l_start_const={0.0397,0.9446,0.0157},
      T_l_profile=true,
      T_v_profile=true,
      redeclare package MediumVapour = 
          ApplicationsThermalSeparation.Media.IdealGasMixtures.H2O_CO2_forPI,
      redeclare package MediumLiquid = 
          ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS_WaterProp,
      redeclare model Reaction = 
          ApplicationsThermalSeparation.Obsolete.CO2_Siemens (
          moleBalance=true,
          a_R={0},
          b_R={-50e3}),
      n_elements=30,
      eps_liq_start=0.05*ones(30),
      redeclare model PressureLoss = 
          ThermalSeparation.PressureLoss.StructuredPackedColumn.Stichlmair (
                                                                  rho_const=false,
            rho_v_nom=0.7),
      redeclare record Geometry = 
          ThermalSeparation.Geometry.StructuredPackedColumn.Mellapak250Y (
          H=15,
          zeta=10,
          d=0.15),
      p_v_start_inlet=105000,
      p_v_start_outlet=100001,
      T_vap_start_bottom=373.15,
      T_vap_start_top=371.15,
      T_liq_start_bottom=373.15,
      T_liq_start_top=371.15,
      T_vapour_start=373.15,
      T_liquid_start=373.15,
      redeclare model FilmModel = 
          ThermalSeparation.FilmModel.StructuredPackedColumn.MS (
            redeclare model HeatTransferVapour = 
              ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Vapour.ChiltonColburn,
            redeclare model HeatTransferLiquid = 
              ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Liquid.ChiltonColburn))
                                             annotation (Placement(transformation(
            extent={{134,98},{182,144}},
                                       rotation=0)));

    Modelica.Blocks.Sources.Ramp ramp2(
      startTime=0,
      duration=1,
      height=0,
      offset=0.44) annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=0,
          origin={78,10})));
    ApplicationsThermalSeparation.Components.HeatExchanger.LiquidTube2
      liquidTube(
      height_in=0,
      zeta=0.02,
      redeclare package MediumLiquid = 
          ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS_WaterProp)
                                                                           annotation (Placement(transformation(
          extent={{-22.5,-20.5},{22.5,20.5}},
          rotation=90,
          origin={200.5,-18.5})));
  Real k=1e-1;
  Real omega = min(1,1 + tanh(k*(time   - 500)));

    SourcesSinks.SinkLiquid sinkLiquid(
      redeclare package Medium = 
          ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS) 
      annotation (Placement(transformation(extent={{-30,-60},{-10,-40}})));
  ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid1(
      x={0,1,0},
      redeclare package MediumLiquid = 
          ApplicationsThermalSeparation.Media.WaterBasedLiquid.BCBO_H2O_AAS_WaterProp,
      T=318.15)                 annotation (Placement(transformation(extent={{82,162},
              {102,182}},   rotation=0)));

    inner SystemTS systemTS 
      annotation (Placement(transformation(extent={{86,86},{106,106}})));
  equation

    connect(reboiler.heatPort,prescribedHeatFlow1. port) annotation (Line(
        points={{227,38.62},{234.5,38.62},{234.5,38},{242,38}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(Desorber.upStreamOut,sinkGas1. gasPortIn) annotation (Line(
        points={{143.6,141.7},{131.8,141.7},{131.8,161},{122,161}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(reboiler.liquidOut,sump. portInRecirc) annotation (Line(
        points={{216.5,54.26},{216.5,82.2},{181.3,82.2}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(sump.portIn,Desorber. downStreamOut) annotation (Line(
        points={{173.2,87.8},{173.2,94.9},{172.4,94.9},{172.4,100.3}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(Desorber.heatPort,ambientHeatSink2. heatPort1) annotation (Line(
        points={{172.88,123.3},{188,126},{191.6,126.8}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(sump.portOut,idealPump1. liquidIn) annotation (Line(
        points={{172.84,73.2},{132,73.2},{132,-44.4}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(ramp3.y, prescribedHeatFlow1.Q_flow) annotation (Line(
        points={{277,38},{262,38}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sump.portOut, liquidTube.liquidIn) annotation (Line(
        points={{172.84,73.2},{182,-14},{180.41,-18.5}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(liquidTube.liquidOut, reboiler.liquidIn) annotation (Line(
        points={{220.59,-18.5},{215,-18.5},{215,27.4}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));

    connect(reboiler.vapourOut, sump.gasPortIn) annotation (Line(
        points={{212.9,54.26},{212.9,66},{158,66},{158,85},{169.24,85}},
        color={160,160,0},
        thickness=1,
        smooth=Smooth.None));
    connect(sump.gasPortOut, Desorber.upStreamIn) annotation (Line(
        points={{170.5,87.8},{158.25,87.8},{158.25,100.3},{143.6,100.3}},
        color={160,160,0},
        thickness=1,
        smooth=Smooth.None));
    connect(sinkLiquid.liquidPortIn, idealPump1.liquidOut) annotation (Line(
        points={{-20,-42.8},{52,-42.8},{52,-63.36},{125.04,-63.36}},
        color={0,80,160},
        thickness=1,
        smooth=Smooth.None));
    connect(ramp2.y, idealPump1.VdotInput) annotation (Line(
        points={{86.8,10},{86.8,-17},{126.72,-17},{126.72,-46.56}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sourceLiquid1.liquidPortOut, Desorber.downStreamIn) annotation (Line(
        points={{92,163.4},{92,154},{172.4,154},{172.4,141.7}},
        color={0,80,160},
        thickness=1,
        smooth=Smooth.None));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-110},
              {200,200}}),       graphics), Icon(coordinateSystem(
            preserveAspectRatio=true, extent={{-100,-110},{200,200}})),
      experiment(StopTime=18000, NumberOfIntervals=5500),
      experimentSetupOutput(equdistant=false));
  end Control3_2_tank_torricelli;

end Tests;
