within ThermalSeparation.Components.LiquidVolumes;
package Tests
  extends Icons.Library.Red;
  model SumpStandalone

  replaceable package MediumLiquid =
  ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O      constrainedby
      ThermalSeparation.Media.BaseMediumLiquid                                                          annotation(choicesAllMatching);
     outer ThermalSeparation.SystemTS systemTS;
   parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

   MediumLiquid.BaseProperties mediumLiquid(T0=T_ref,p=p, T=T, x=x_l);
      MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref,p=1e5, T=T_in, x=x_in);
         MediumLiquid.BaseProperties mediumLiquidInRecirc(T0=T_ref,p=p_in, T=T_in_recirc, x=x_in_recirc);

      parameter Boolean inertLiquid[nSL] = fill(false,nSL);
    parameter Integer nS=3
      "number of species which are equal in vapour and liquid phase";
    final parameter Integer nL=nSL-nS
      "number of additional substances which are only in liquid phase";

    final parameter Integer nSL = MediumLiquid.nSubstance;

  /*** Medium properties ***/
    SI.Density rho_l_in = mediumLiquidIn.d;
    SI.Concentration c_l_in[nSL];// = portIn.c;
    ThermalSeparation.Units.MolarEnthalpy h_l_in= mediumLiquidIn.h;
    SI.MolarMass MM_l_in = mediumLiquidIn.MM;

    SI.Density rho_l_in_recirc = mediumLiquidInRecirc.d;
    SI.Concentration c_l_in_recirc[nSL]; //= portInRecirc.c;
    ThermalSeparation.Units.MolarEnthalpy h_l_in_recirc= mediumLiquidInRecirc.h;
    SI.MolarMass MM_l_in_recirc = mediumLiquidInRecirc.MM;

    SI.Density rho_l = mediumLiquid.d;
    SI.Concentration c_l[nSL];
     ThermalSeparation.Units.MolarEnthalpy h_l= mediumLiquid.h;
     ThermalSeparation.Units.MolarEnthalpy u_l = mediumLiquid.u;
     SI.MolarMass MM_l = mediumLiquid.MM;
    SI.MoleFraction x_l[nSL];

   SI.Temperature T;

   SI.VolumeFlowRate Vdot_l_out(start=1);
   SI.VolumeFlowRate Vdot_l_in;// = portIn.Recirc.Vdot;
   SI.VolumeFlowRate Vdot_l_in_recirc; //= portInRecirc.Vdot;

    SI.Pressure p(start=2e5);

    /*** geometry data ***/
      final parameter SI.Area A= Modelica.Constants.pi/4* d_volume^2;
    SI.Height level(stateSelect=StateSelect.always);

    parameter SI.Length height_in = -0.15;
    parameter SI.Length length_out = 0.15;
    parameter SI.Diameter d_volume = 0.025;
    parameter Real zeta=2;

    parameter SI.Height inletHeight = 0;
    parameter SI.Pressure p_ambient = 1e5;

  /*** Initialization ***/
  parameter SI.Temperature T_start=300;
  parameter Real x_start[nSL] = {0.0226,1-0.0226-0.0313,0.0313};
  parameter SI.Height level_start= 0.1;

  //   Interfaces.LiquidPortOut portOut(nS=nSL)
  //     annotation (Placement(transformation(extent={{-34,-78},{-14,-58}})));
  //   Interfaces.LiquidPortIn portInRecirc(nS=nSL)
  //     annotation (Placement(transformation(extent={{60,12},{80,32}})));
  //
  //     Interfaces.LiquidPortIn portIn(nS=nSL)
  //     annotation (Placement(transformation(extent={{-30,68},{-10,88}})));
  //
  //   Modelica.Blocks.Interfaces.RealOutput y=level
  //     annotation (Placement(transformation(extent={{74,58},{94,78}})));

  /*** monitoring ***/
  SI.MassFlowRate mdot_in = Vdot_l_in*1000;
  SI.MassFlowRate mdot_in2=Vdot_l_in_recirc*1000;
  SI.MassFlowRate mdot_out=Vdot_l_out*1000;
  Real sum_x = sum(x_l);

  /*** parameters to be set if model is used standalone ***/
  Real T_in = 95+273;
  Real T_in_recirc = 96+273;
  Real x_in[nSL]={0.05,0.9,0.05};
  Real x_in_recirc[nSL]={0.05,0.9,0.05};
  Real p_in=1e5;
  Real p_out = 1.3e5;

  equation
    Vdot_l_in = 0.0003;
    Vdot_l_in_recirc=0.0003;
    p_out=p;
    c_l_in = x_in*rho_l_in/MM_l_in;
     c_l_in_recirc = x_in_recirc*rho_l_in_recirc/MM_l_in_recirc;
  //   portOut.c = c_l;
  //   portOut.x = x_l;
  //   portOut.p = p;
  //   portOut.T = T;
  //   portOut.Vdot = - Vdot_l_out;

    /*** correlation between mole fraction x and concentration c ***/
     for i in 1:nSL loop
       x_l[i] = c_l[i] * MM_l /rho_l;
     end for;

       /***energy balance ***/
     A* der(sum(c_l[:])*u_l*level)  =   Vdot_l_in*sum(c_l_in[:])*h_l_in + Vdot_l_in_recirc*sum(c_l_in_recirc[:])*h_l_in_recirc - Vdot_l_out*sum(c_l[:])*h_l;
       /*** mass balance ***/
      for i in 1:nSL loop
              A* der(c_l[i]*level) = Vdot_l_in*c_l_in[i] + Vdot_l_in_recirc*c_l_in_recirc[i] - Vdot_l_out*c_l[i];
      end for;
         A* der(rho_l/MM_l*level) = Vdot_l_in*rho_l_in/MM_l_in + Vdot_l_in_recirc*rho_l_in_recirc/MM_l_in_recirc - Vdot_l_out*rho_l/MM_l;
     // sum(x_l)=1;

  /*** fill level ***/
  //A*der(level) = Vdot_l_in + Vdot_l_in_recirc - Vdot_l_out;

  p=level*9.8*rho_l+p_in- zeta* rho_l/2*((Vdot_l_in+Vdot_l_in_recirc)/A)^2;

  initial equation
     T=T_start;
     x_l=x_start;

   // level=level_start;

   // p=1.015e5;

    annotation (Diagram(graphics));
  end SumpStandalone;

  model TestSumpStandalone

    SumpStandalone sumpStandalone
      annotation (Placement(transformation(extent={{-54,8},{-34,28}})));
    inner SystemTS systemTS
      annotation (Placement(transformation(extent={{16,40},{36,60}})));
    annotation (Commands(file="../../ttt.mos" "ttt"));
  end TestSumpStandalone;

  model TestTank

    ThermalSeparation.Components.LiquidVolumes.Tank2Inlets tank
      annotation (Placement(transformation(extent={{-4,-10},{16,10}})));
    SourcesSinks.SourceLiquid sourceLiquid(redeclare package MediumLiquid =
          ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O, x={0,1})
      annotation (Placement(transformation(extent={{-40,16},{-20,36}})));
    SourcesSinks.SourceLiquid sourceLiquid1(
      redeclare package MediumLiquid =
          ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
      x={0,1},
      fixed_pressure=false)
      annotation (Placement(transformation(extent={{14,26},{34,46}})));
    SourcesSinks.SinkLiquid sinkLiquid(
      redeclare package Medium =
          ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
      use_p=false,
      use_Vdot=true)
      annotation (Placement(transformation(extent={{-8,-50},{12,-30}})));
    inner SystemTS systemTS
      annotation (Placement(transformation(extent={{-44,66},{-24,86}})));
    Modelica.Blocks.Sources.Ramp coolingWater1(
      duration=10,
      height=0,
      offset=2e-4,
      startTime=0)
      annotation (Placement(transformation(extent={{-48,-50},{-32,-34}})));
  equation
    connect(sinkLiquid.liquidPortIn, tank.portOut) annotation (Line(
        points={{2,-32.8},{4,-32.8},{4,-8.2},{6.2,-8.2}},
        color={0,80,160},
        thickness=1,
        smooth=Smooth.None));
    connect(sourceLiquid.liquidPortOut, tank.portIn) annotation (Line(
        points={{-30,17.4},{-14,17.4},{-14,8.4},{2.4,8.4}},
        color={0,80,160},
        thickness=1,
        smooth=Smooth.None));
    connect(sourceLiquid1.liquidPortOut, tank.portInRecirc) annotation (Line(
        points={{24,27.4},{18,27.4},{18,8.4},{9.4,8.4}},
        color={0,80,160},
        thickness=1,
        smooth=Smooth.None));
    connect(coolingWater1.y, sinkLiquid.Vdot_set) annotation (Line(
        points={{-31.2,-42},{-18,-42},{-18,-36.7},{-5.1,-36.7}},
        color={0,0,127},
        smooth=Smooth.None));
    annotation (Diagram(graphics));
  end TestTank;

  model TestSump

    Sump tank(level_start=0.1)
      annotation (Placement(transformation(extent={{-4,-10},{16,10}})));
    SourcesSinks.SourceLiquid sourceLiquid(redeclare package MediumLiquid =
          ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O, x={0,1})
      annotation (Placement(transformation(extent={{-22,42},{-2,62}})));
    SourcesSinks.SourceLiquid sourceLiquid1(
      redeclare package MediumLiquid =
          ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
      x={0,1},
      fixed_pressure=false)
      annotation (Placement(transformation(extent={{14,26},{34,46}})));
    SourcesSinks.SinkLiquid sinkLiquid(
      redeclare package Medium =
          ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O,
      use_Vdot=true,
      use_p=true,
      p=100000)
      annotation (Placement(transformation(extent={{30,-62},{50,-42}})));
    inner SystemTS systemTS
      annotation (Placement(transformation(extent={{-44,66},{-24,86}})));
    SourcesSinks.SourceGas sourceGas(redeclare package Medium =
          ThermalSeparation.Media.IdealGasMixtures.H2O_CO2, use_Flow=false)
      annotation (Placement(transformation(extent={{-38,-12},{-18,8}})));
    SourcesSinks.SinkGas sinkGas(redeclare package Medium =
          ThermalSeparation.Media.IdealGasMixtures.H2O_CO2, use_p=false)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-30,30})));
    Modelica.Blocks.Sources.Ramp ramp(
      height=0,
      offset=1,
      duration=0)
      annotation (Placement(transformation(extent={{0,-78},{20,-58}})));
  equation
    connect(sourceLiquid.liquidPortOut, tank.portIn) annotation (Line(
        points={{-0.6,52},{6,52},{6,9.6}},
        color={0,80,160},
        thickness=1,
        smooth=Smooth.None));
    connect(sourceLiquid1.liquidPortOut, tank.portInRecirc) annotation (Line(
        points={{35.4,36},{42,36},{42,2},{14.6,2}},
        color={0,80,160},
        thickness=1,
        smooth=Smooth.None));
    connect(sinkGas.gasPortIn, tank.gasPortOut) annotation (Line(
        points={{-20.4,30},{-6,30},{-6,7.8},{-2.4,7.8}},
        color={160,160,0},
        thickness=1,
        smooth=Smooth.None));
    connect(sourceGas.gasPortOut, tank.gasPortIn) annotation (Line(
        points={{-16.6,-2},{-14.1,-2},{-14.1,3.6},{-2.4,3.6}},
        color={160,160,0},
        thickness=1,
        smooth=Smooth.None));
    connect(sinkLiquid.liquidPortIn, tank.portOut) annotation (Line(
        points={{30.4,-52},{6,-52},{6,-9.6}},
        color={153,217,234},
        thickness=1));
    connect(ramp.y, sinkLiquid.Vdot_set) annotation (Line(points={{21,-68},{38,
            -68},{38,-70},{56,-70},{56,-42},{48.5,-42},{48.5,-43.5}}, color={0,
            0,127}));
  end TestSump;
end Tests;
