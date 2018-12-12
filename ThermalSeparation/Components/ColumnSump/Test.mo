within ThermalSeparation.Components.ColumnSump;
package Test
  model Test_SumpSteam

    Sump_Steam_1_2 sump_Steam(
      redeclare package MediumVapour = 
          ThermalSeparation.Media.C2H5OH_Water_Vap,
      redeclare package MediumLiquid = 
          ThermalSeparation.Media.C2H5OH_Water_Liq,
      vdot_steam_in=1,
      redeclare model HeatTransferModel = 
          ThermalSeparation.Components.ColumnSump.HeatTransferModel.dT_m_2) 
      annotation (Placement(transformation(extent={{-58,-20},{-6,36}})));
    ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiq_x2_HEX(
      x={0.1,0.9},
      redeclare package MediumLiquid = 
          ThermalSeparation.Media.C2H5OH_Water_Liq,
      nS=2,
      T=343.15,
      Flow=0.001,
      p=120000) 
      annotation (Placement(transformation(extent={{-30,46},{-10,66}})));

    Components.SourcesSinks.SinkGas sinkGas1(          nS=2, p=100000) 
      annotation (Placement(transformation(extent={{-62,48},{-42,68}})));
    Components.SourcesSinks.SinkLiquid sinkLiquid1(          nS=2, use_p=
          false) 
      annotation (Placement(transformation(extent={{-45,-44},{-25,-24}})));
      inner SystemTS systemTS(                  T_ref=273.15) 
        annotation (Placement(transformation(extent={{24,72},{44,92}})));
  equation
    connect(sourceLiq_x2_HEX.liquidPortOut, sump_Steam.liquidIn) annotation (
        Line(
        points={{-20,47.4},{-20,40},{-27.32,40},{-27.32,27.6}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(sinkGas1.gasPortIn, sump_Steam.vapourOut) annotation (Line(
        points={{-52,51},{-52,40},{-42.92,40},{-42.92,27.6}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(sump_Steam.liquidOut, sinkLiquid1.liquidPortIn) annotation (Line(
        points={{-35.12,-7.68},{-35.12,-20.12},{-35,-20.12},{-35,-26.8}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    annotation (Diagram(graphics));
  end Test_SumpSteam;

  model Test_SumpSteam_2_4

    Sump_Steam_2_4 sump_Steam(redeclare model HeatTransferModel = 
          ThermalSeparation.Components.ColumnSump.HeatTransferModel.dT_m_2) 
      annotation (Placement(transformation(extent={{-58,-20},{-6,36}})));
    Components.SourcesSinks.SinkGas sinkGas1(          nS=2, p=100000) 
      annotation (Placement(transformation(extent={{-58,44},{-38,64}})));
    Components.SourcesSinks.SinkLiquid sinkLiquid1(          nS=2, use_p=
          false) 
      annotation (Placement(transformation(extent={{-45,-38},{-25,-18}})));
    ThermalSeparation.Kai.SourcesSinksxx.SinkGas sinkGas2(
                                  nS=1) 
      annotation (Placement(transformation(extent={{-72,18},{-52,38}})));
    ThermalSeparation.Kai.SourcesSinksxx.SinkLiquid_x sinkLiquid2( use_p=
          false, nS=1) 
      annotation (Placement(transformation(extent={{-72,-22},{-52,-2}})));
    ThermalSeparation.Kai.SourcesSinksxx.SourceGas_Vdot_IF97
      sourceGas_Vdot_IF97_1(
      nS=1,
      T=393.15,
      p=198000) annotation (Placement(transformation(extent={{-8,-8},{12,12}})));
    Modelica.Blocks.Sources.RealExpression Vdot(y=0.5) 
      annotation (Placement(transformation(extent={{12,-28},{-8,-9}})));
      inner SystemTS systemTS(                  T_ref=273.15) 
        annotation (Placement(transformation(extent={{14,62},{34,82}})));
    ThermalSeparation.Components.SourcesSinks.SourceLiquid sourceLiquid_x(
      redeclare package MediumLiquid = 
          ThermalSeparation.Media.C2H5OH_Water_Liq,
      nS=2,
      x={0.1,0.9},
      T=348.15,
      Flow=0.001) 
      annotation (Placement(transformation(extent={{-28,44},{-8,64}})));

  equation
    connect(sinkGas1.gasPortIn, sump_Steam.vapourOut) annotation (Line(
        points={{-48,47},{-48,38},{-42.92,38},{-42.92,27.6}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(sump_Steam.liquidOut, sinkLiquid1.liquidPortIn) annotation (Line(
        points={{-35.12,-7.68},{-35.12,-20.12},{-35,-20.12},{-35,-20.8}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(sinkGas2.gasPortIn, sump_Steam.gasPortOut) annotation (Line(
        points={{-62.2,18.6},{-62.2,13.88},{-50.46,13.88}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(sinkLiquid2.liquidPortIn, sump_Steam.liquidPortOut) annotation (
        Line(
        points={{-62,-4},{-62,0.44},{-50.46,0.44}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(sourceGas_Vdot_IF97_1.gasPortOut, sump_Steam.gasPortIn) annotation (
       Line(
        points={{2,10},{2,13.88},{-19.78,13.88}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(Vdot.y, sourceGas_Vdot_IF97_1.Vdot) annotation (Line(
        points={{-9,-18.5},{-13.5,-18.5},{-13.5,3.6},{-7.4,3.6}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sourceLiquid_x.liquidPortOut, sump_Steam.liquidIn) annotation (Line(
        points={{-18,45.4},{-18,38},{-27.32,38},{-27.32,27.6}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    annotation (Diagram(graphics));
  end Test_SumpSteam_2_4;
end Test;
