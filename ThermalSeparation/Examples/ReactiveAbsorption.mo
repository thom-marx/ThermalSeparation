within ThermalSeparation.Examples;
model ReactiveAbsorption
  "absorption with reaction of SO2 in liquid phase, spray absorber"
      //Lauftest: 1.2.
  //Lauftest: 21.12.

ThermalSeparation.Components.SourcesSinks.SourceLiquid
                                              sourceLiquid(
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O_H_HSO3,
    x={0,0,0,0,0,0,1,0,0},
    Flow=0.03,
    T=323.15)                          annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
                                   rotation=270,
        origin={12,56})));

ThermalSeparation.Components.SourcesSinks.SinkLiquid
                                            sinkLiquid(
                                                     redeclare package Medium =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O_H_HSO3)
                      annotation (Placement(transformation(extent={{-10,-10},{
            10,10}},
                   rotation=270,
        origin={20,-48})));

ThermalSeparation.Components.SourcesSinks.SinkGas
                                         sinkGas(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_SO2_HCl_HF_N2,
    p=149910)   annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-30,56})));

ThermalSeparation.Components.Columns.SprayColumn column(
    n_elements=8,
    hasLiquidFeed=false,
    hasVapourFeed=false,
    redeclare record Geometry =
        ThermalSeparation.Geometry.StructuredPackedColumn.Geometry (d=4.8, H=10),
    redeclare model HeatTransferWall = ThermalSeparation.Wall.Adiabatic,
    nS=7,
    inertVapour={false,false,false,false,false,false,false},
    mapping={{1,7},{2,2},{3,3},{4,4},{5,5},{6,6},{7,1}},
    redeclare package MediumVapour =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_SO2_HCl_HF_N2,
    redeclare package MediumLiquid =
        ThermalSeparation.Media.WaterBasedLiquid.N2_O2_CO2_SO2_HCl_HF_H2O_H_HSO3,
    inertLiquid={false,false,false,false,false,false,false,true,true},
    x_l_start_const={1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1 - 8e-7,1e-7,1e-7},
    redeclare model Reaction =
        ThermalSeparation.Reaction.ReactionEquilibrium.TrueEquilibrium,
    redeclare model InitOption =
        ThermalSeparation.Components.Columns.BaseClasses.Initialization.Special_x_l,
    x_v_start_const={0,0,0.0496,0.0494,0.004859,0.004857,0.7838},
    c_v(start=fill({1,1,1,1,1,1,60}, 8)),
    c_l(start=fill({500,500,500,500,500,500,50000,1,1}, 8)),
    redeclare model BalanceEquations =
        ThermalSeparation.BalanceEquations.SprayColumn.NonEquilibrium.TwoPhaseVarState (
         redeclare model FilmModel = ThermalSeparation.FilmModel.SprayColumn.MS (
             redeclare replaceable model StateSelection =
                ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.None)),
    p_v_start_inlet=150500,
    p_v_start_outlet=150100,
    T_vapour_start(displayUnit="degC") = 323.15,
    T_liquid_start(displayUnit="degC") = 323.15)
                                     annotation (Placement(transformation(extent={{-34,-14},
            {14,32}}, rotation=0)));

  ThermalSeparation.Components.SourcesSinks.SourceGas
                                             sourceGas_Vdot(
    redeclare package Medium =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_SO2_HCl_HF_N2,
    x={0.05,0.05,0.05,0.05,0.005,0.005,0.79},
    T=417.15)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=90,
        origin={-30,-48})));

  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    offset=120,
    height=-60,
    startTime=200)
              annotation (Placement(transformation(extent={{-68,-64},{-48,-44}},
          rotation=0)));
  ThermalSeparation.Components.SourcesSinks.AmbientHeatSink ambientHeatSink
    annotation (Placement(transformation(extent={{28,-6},{48,14}},rotation=0)));

  inner SystemTS systemTS
    annotation (Placement(transformation(extent={{-64,2},{-44,22}})));
  Modelica.Blocks.Sources.Constant const(k=0.5)
    annotation (Placement(transformation(extent={{-14,-58},{6,-38}})));
equation
  connect(ramp.y, sourceGas_Vdot.Flow_in) annotation (Line(
      points={{-47,-54},{-44,-54},{-44,-58},{-34,-58}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sinkLiquid.liquidPortIn, column.downStreamOut) annotation (Line(
      points={{20,-38.4},{20,-26},{6.8,-26},{6.8,-11.7}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(sourceGas_Vdot.gasPortOut, column.upStreamIn) annotation (Line(
      points={{-30,-36.6},{-30,-26},{-26.8,-26},{-26.8,-11.7}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  connect(ambientHeatSink.heatPort1, column.heatPort) annotation (Line(
      points={{26.4,4},{20,4},{20,9},{9.68,9}},
      color={188,51,69},
      smooth=Smooth.None));
  connect(column.downStreamIn, sourceLiquid.liquidPortOut) annotation (Line(
      points={{6.8,29.7},{6.8,36},{12,36},{12,44.6}},
      color={153,217,234},
      thickness=1,
      smooth=Smooth.None));
  connect(sinkGas.gasPortIn, column.upStreamOut) annotation (Line(
      points={{-30,46.4},{-30,36},{-26.8,36},{-26.8,29.7}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                    graphics),
    experiment(StopTime=100, Algorithm="Dassl"),
    experimentSetupOutput,
    Documentation(info="<html>
<p><ul>
<li> spray column</li>
<li> 7 components in vapour phase </li>
<li> 9 components in liquid phase (2 inert components)</li>
<li> absorption process</li>
<li> reaction in liquid phase</li>
</ul></p>
</html>"));
end ReactiveAbsorption;
