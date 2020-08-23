within ThermalSeparation.Components;
package HeatExchanger
extends ThermalSeparation.Icons.Library.Red;

  model LiquidCooler "cools fluid so specified temperature, no pressure loss"
  extends ThermalSeparation.Icons.Color.HeatExchanger2;
    //import CCS.ThermalSeparation;

  Modelica.SIunits.HeatFlowRate Q;
  parameter Modelica.SIunits.Temperature T_set "temperature at outlet";

  outer ThermalSeparation.SystemTS systemTS;

  /*** Medium properties ***/

  MediumLiquid.BaseProperties medium(p=coldLiquidOut.p, T=T_set, x=coldLiquidOut.x_outflow,h=h);

  replaceable package MediumLiquid =
       ThermalSeparation.Media.WaterBasedLiquid.N2_H2O
    constrainedby ThermalSeparation.Media.BaseMediumLiquid "medium to be used" annotation(choicesAllMatching);

    ThermalSeparation.Interfaces.LiquidPortIn hotLiquidIn(redeclare package Medium =
                 MediumLiquid) annotation (Placement(transformation(extent={{-80,
              54},{-60,74}}), iconTransformation(extent={{-120,14},{-80,54}})));
    ThermalSeparation.Interfaces.LiquidPortOut coldLiquidOut(
      redeclare package Medium = MediumLiquid,
      h_outflow=h,
      x_outflow=x_in,
      p=hotLiquidIn.p) annotation (Placement(transformation(extent={{62,52},{82,
              72}}), iconTransformation(extent={{-120,-52},{-80,-12}})));
    ThermalSeparation.Interfaces.HeatPort heatPort(Qdot=-Q) annotation (Placement(
          transformation(extent={{82,0},{102,20}}, rotation=0),
          iconTransformation(extent={{62,-20},{102,20}})));
  Modelica.SIunits.MoleFraction x_in[MediumLiquid.nSubstance];
  Real h_in;
  Real h;//=medium.h;
  equation
          hotLiquidIn.x_outflow= inStream(coldLiquidOut.x_outflow);
           hotLiquidIn.h_outflow= inStream(coldLiquidOut.h_outflow);
  x_in = inStream(hotLiquidIn.x_outflow);
  h_in = inStream(hotLiquidIn.h_outflow);

  /* energy balance */
  Q =hotLiquidIn.Ndot*(coldLiquidOut.h_outflow - h_in);// (-coldLiquidOut.Vdot)*mediumLiquid.d*calcSpecificEnthalpy.cp_tot*(hotLiquidIn.T-T_set);

  /* mass balance */
  hotLiquidIn.Ndot + coldLiquidOut.Ndot=0;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                     graphics), Diagram(coordinateSystem(preserveAspectRatio=
              false, extent={{-100,-100},{100,100}}),
                                        graphics));
  end LiquidCooler;

  model CcounterFlowHeatExchanger "ideal - no pressure loss"
       extends ThermalSeparation.Icons.Color.HeatExchanger;
    import SI = Modelica.SIunits;

  replaceable package MediumLiquidCold =
        ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O
    constrainedby ThermalSeparation.Media.BaseMediumLiquid "medium to be used for hot liquid" annotation(choicesAllMatching);
  replaceable package MediumLiquidHot =
        ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O
    constrainedby ThermalSeparation.Media.BaseMediumLiquid "medium to be used for cold liquid" annotation(choicesAllMatching);

  // parameter
      // Geometry
    parameter SI.Mass plateMetalMass = 310 annotation(Dialog(group="Geometry"));
    parameter SI.SpecificHeatCapacity cpMetal = 500 annotation(Dialog(group="Geometry"));
    parameter SI.Volume volumeHotLiquid = 11e-3 annotation(Dialog(group="Geometry"));
    parameter SI.Volume volumeColdLiquid = 11e-3 annotation(Dialog(group="Geometry"));

    parameter SI.Area surfaceArea = 10.4 annotation(Dialog(group="Geometry"));
    parameter Boolean alpha_input=false;
     SI.CoefficientOfHeatTransfer alphaHot_nom;
     SI.CoefficientOfHeatTransfer alphaCold_nom;
     SI.CoefficientOfHeatTransfer alphaHot_nom_param = 3600 annotation(Dialog(group="Geometry"));
     SI.CoefficientOfHeatTransfer alphaCold_nom_param = 3600 annotation(Dialog(group="Geometry"));

      // nominal operating point
    parameter SI.Pressure dpHotLiquid = 2e5 annotation(Dialog(group="nominal operating point"));
    parameter SI.Pressure dpColdLiquid = 2e5 annotation(Dialog(group="nominal operating point"));
    parameter SI.MassFlowRate mFlowHotLiquid_nom = 0.3 annotation(Dialog(group="nominal operating point"));
    parameter SI.MassFlowRate mFlowColdLiquid_nom = 0.3 annotation(Dialog(group="nominal operating point"));
    parameter SI.Temperature ThotLiquidOut_nom = ThotLiquid_start annotation(Dialog(group="nominal operating point"));
    parameter SI.Temperature TcoldLiquidOut_nom = TcoldLiquid_start annotation(Dialog(group="nominal operating point"));

    parameter SI.AbsolutePressure pLiquid_nom = 4e5 annotation(Dialog(group="nominal operating point"));

      // initial
    parameter SI.Temperature ThotLiquidIn_start = 104 + 273.15 annotation(Dialog(tab="Initialization"));
    parameter SI.Temperature TcoldLiquidIn_start = 51 + 273.15 annotation(Dialog(tab="Initialization"));
    parameter SI.Temperature ThotLiquid_start = 56 + 273.15 annotation(Dialog(tab="Initialization"));
    parameter SI.Temperature TcoldLiquid_start = 98 + 273.15 annotation(Dialog(tab="Initialization"));

      // Medium parameter
    final parameter Integer nShot = 2 "number of species which are equal in vapour and liquid phase";
    parameter Integer nScold = 2 "number of species which are equal in vapour and liquid phase";

   // final parameter

    final parameter Integer nLhot = nSLhot - nShot "number of additional substances which are only in liquid phase";
    final parameter Integer nLcold = nSLcold - nScold "number of additional substances which are only in liquid phase";
    final parameter Integer nSLhot= MediumLiquidHot.nSubstance;
    final parameter Integer nSLcold = MediumLiquidCold.nSubstance;

    final parameter SI.Density rhoHotLiquid_nom = MediumLiquidHot.DensityWater(p=pLiquid_nom, T=ThotLiquidOut_nom);
    final parameter SI.Density rhoColdLiquid_nom = MediumLiquidCold.DensityWater(p=pLiquid_nom, T=TcoldLiquidOut_nom);

    final parameter SI.Length dummyDiameter = 0.0254 " diameter is abitrary"; // since dp and volume is adjusted using tube zeta resp. length
    final parameter SI.Length hotLiquidTubeLength = volumeHotLiquid/(0.25*Modelica.Constants.pi*dummyDiameter^2);
    final parameter SI.Length coldLiquidTubeLength = volumeColdLiquid/(0.25*Modelica.Constants.pi*dummyDiameter^2);
    final parameter SI.Area CrossSectionArea = Modelica.Constants.pi*dummyDiameter^2/4;
    final parameter Real zetaHotLiquid = (dpHotLiquid*CrossSectionArea^2*2*dummyDiameter)/((mFlowHotLiquid_nom/rhoHotLiquid_nom)^2*hotLiquidTubeLength*rhoHotLiquid_nom);
    final parameter Real zetaHotLiquid_haiko = Modelica.Constants.pi^2*dummyDiameter^5/(8*hotLiquidTubeLength)*(dpHotLiquid*rhoHotLiquid_nom/mFlowHotLiquid_nom^2);
    final parameter Real zetaColdLiquid = Modelica.Constants.pi^2*dummyDiameter^5/(8*coldLiquidTubeLength)*(dpColdLiquid*rhoColdLiquid_nom/mFlowColdLiquid_nom^2);
  //parameter SI.Area areaHT = 1 "heat transfer area";

    parameter SI.Concentration c_l_start_hot[nSLhot]={5,40000} "initial value for concentration of hot liquid" annotation(Dialog(tab="Initialization"));
    parameter SI.Concentration c_l_start_cold[nSLcold]={5,40000} "initial value for concentration of cold liquid" annotation(Dialog(tab="Initialization"));

  //variables
  SI.TemperatureDifference logDtHot;
  SI.TemperatureDifference logDtCold;

  SI.Temperature TwallHot(start=(ThotLiquidIn_start+TcoldLiquid_start)/2);
  SI.Temperature TwallCold(start=(TcoldLiquidIn_start+ThotLiquid_start)/2);

  SI.HeatFlowRate QFlowFromColdLiqidInToColdWall;
  SI.HeatFlowRate QFlowFromColdLiqidOutToHotWall;
  SI.HeatFlowRate QFlowFromHotLiqidInToHotWall;
  SI.HeatFlowRate QFlowFromHotLiqidOutToColdWall;
  //Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatFlowHotLiquid;
  //Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatFlowColdLiquid;
  SI.CoefficientOfHeatTransfer alphaCold;
  SI.CoefficientOfHeatTransfer alphaHot;
  /*** state selection ***/
  parameter Boolean ss_hot_c2=true "true if c_l[2] of the hot liquid is to be a state";
  parameter Boolean ss_cold_c2=true "true if c_l[2] of the cold liquid is to be a state";

  // components

    ThermalSeparation.Components.HeatExchanger.LiquidTube hotLiquidTube(
      redeclare replaceable package MediumLiquid = MediumLiquidHot,
      d=dummyDiameter,
      height_in=0,
      zeta=zetaHotLiquid,
      length_out=hotLiquidTubeLength,
      d_volume=dummyDiameter,
      T_start=ThotLiquid_start,
      useX_start=false,
      ss_c2=ss_hot_c2,
      c_l_start=c_l_start_hot)
      annotation (Placement(transformation(extent={{-70,-20},{-38,20}})));
    ThermalSeparation.Components.HeatExchanger.LiquidTube coldLiquidTube(
      redeclare replaceable package MediumLiquid = MediumLiquidCold,
      d=dummyDiameter,
      height_in=0,
      zeta=zetaColdLiquid,
      length_out=coldLiquidTubeLength,
      d_volume=dummyDiameter,
      T_start=TcoldLiquid_start,
      useX_start=false,
      ss_c2=ss_cold_c2,
      c_l_start=c_l_start_cold)  annotation (Placement(transformation(
          extent={{-16,-20},{16,20}},
          rotation=180,
          origin={38,0})));
    ThermalSeparation.Interfaces.LiquidPortIn hotLiquidIn(redeclare package
        Medium = MediumLiquidHot) annotation (Placement(transformation(extent={{
              60,60},{80,80}}), iconTransformation(extent={{80,16},{120,56}})));
    ThermalSeparation.Interfaces.LiquidPortOut hotLiquidOut(redeclare package
        Medium = MediumLiquidHot) annotation (Placement(transformation(extent={{-100,
              -100},{-80,-80}}), iconTransformation(extent={{-120,-52},{-80,-12}})));
    ThermalSeparation.Interfaces.LiquidPortIn coldLiquidIn(redeclare package
        Medium = MediumLiquidCold) annotation (Placement(transformation(extent={{
              -100,60},{-80,80}}), iconTransformation(extent={{-120,14},{-80,54}})));
    ThermalSeparation.Interfaces.LiquidPortOut coldLiquidOut(redeclare package
        Medium = MediumLiquidCold) annotation (Placement(transformation(extent={{
              60,-100},{80,-80}}), iconTransformation(extent={{80,-52},{120,-12}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
      heatFlowColdLiquidTube
      annotation (Placement(transformation(extent={{0,-12},{20,8}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
      heatFlowHotLiquidTube
      annotation (Placement(transformation(extent={{-20,-10},{-38,8}})));
  MediumLiquidCold.BaseProperties mediumColdIn(p=coldLiquidIn.p, x=inStream(coldLiquidIn.x_outflow), T=T_c_in,h=h_c_in);
  Real T_c_in;
  Real h_c_in;

  MediumLiquidHot.BaseProperties mediumHotIn(p=hotLiquidIn.p, x=inStream(hotLiquidIn.x_outflow), T=T_h_in,h=h_h_in);
  Real T_h_in;
  Real h_h_in;

  MediumLiquidCold.BaseProperties mediumColdOut(p=coldLiquidOut.p, x=(coldLiquidOut.x_outflow), T=T_c_out,h=h_c_out);
  Real T_c_out;
  Real h_c_out;

  MediumLiquidHot.BaseProperties mediumHotOut(p=hotLiquidOut.p, x=(hotLiquidOut.x_outflow), T=T_h_out,h=h_h_out);
  Real T_h_out;
  Real h_h_out;

  SI.VolumeFlowRate Vdot_c_in = (coldLiquidIn.Ndot)*mediumColdIn.MM/mediumColdIn.d;
  SI.VolumeFlowRate Vdot_h_in = (hotLiquidIn.Ndot)*mediumHotIn.MM/mediumHotIn.d;

  Real omega;
    Modelica.Blocks.Interfaces.RealInput alpha_in if alpha_input annotation (Placement(transformation(extent={{-136,-20},{-96,20}})));

  protected
  Modelica.Blocks.Interfaces.RealInput alpha_internal;

  initial equation
     TwallHot = (ThotLiquidIn_start+TcoldLiquid_start)/2;
     TwallCold = (TcoldLiquidIn_start+ThotLiquid_start)/2;

  equation
    if not alpha_input then
      alpha_internal=3600;
    end if;

    omega = 0.5 + 0.5*tanh(40000*(Vdot_h_in - 3e-4));

    if alpha_input then
       if Vdot_h_in<=4e-4 then
  //      alphaHot_nom=0;
  //      alphaCold_nom=0;
       alphaHot_nom=omega*alpha_internal;
       alphaCold_nom=omega*alpha_internal;
     else
       alphaHot_nom=alpha_internal;
       alphaCold_nom=alpha_internal;
     end if;
    else
     if Vdot_h_in<=4e-4 then
  //      alphaHot_nom=0;
  //      alphaCold_nom=0;
       alphaHot_nom=omega*alphaHot_nom_param;
       alphaCold_nom=omega*alphaCold_nom_param;
     else
       alphaHot_nom=alphaHot_nom_param;
       alphaCold_nom=alphaCold_nom_param;
     end if;
    end if;

      h_c_in = inStream(coldLiquidIn.h_outflow);
      h_h_in = inStream(hotLiquidIn.h_outflow);
      h_c_out = coldLiquidOut.h_outflow;
      h_h_out = hotLiquidOut.h_outflow;

  alphaCold  = alphaCold_nom*(rhoColdLiquid_nom*Vdot_c_in)/mFlowColdLiquid_nom;
  alphaHot = alphaHot_nom*(rhoHotLiquid_nom*Vdot_h_in)/mFlowHotLiquid_nom;
  //  if  (hotLiquidOut.T - TwallCold) < 0.1  or (hotLiquidIn.T - TwallHot) < 0.1 then
  //    logDtHot = ((hotLiquidOut.T - TwallCold) + (hotLiquidIn.T - TwallHot))/2;
  //  else
  //    logDtHot = ((hotLiquidOut.T - TwallCold)-(hotLiquidIn.T - TwallHot))/Modelica.Math.log((hotLiquidOut.T - TwallCold)/(hotLiquidIn.T - TwallHot));
  //  end if;

  //logDtHot = ((hotLiquidOut.T - TwallCold)-(hotLiquidIn.T - TwallHot))/Modelica.Math.log((hotLiquidOut.T - TwallCold)/(hotLiquidIn.T - TwallHot));
  logDtHot = ((T_h_out - TwallCold) + (T_h_in - TwallHot))/2;

  //  if (coldLiquidOut.T - TwallHot) > -0.1 or (coldLiquidIn.T - TwallCold) > -0.1 then
  //      logDtCold =(coldLiquidOut.T - TwallHot)+(coldLiquidIn.T - TwallCold)/2;
  //  else
  //    logDtCold = ((coldLiquidOut.T - TwallHot)-(coldLiquidIn.T - TwallCold))/Modelica.Math.log((coldLiquidOut.T - TwallHot)/(coldLiquidIn.T - TwallCold));
  //  end if;

  //logDtCold = ((coldLiquidOut.T - TwallHot)-(coldLiquidIn.T - TwallCold))/Modelica.Math.log((coldLiquidOut.T - TwallHot)/(coldLiquidIn.T - TwallCold));
  logDtCold =(T_c_out - TwallHot)+(T_c_in - TwallCold)/2;

  heatFlowHotLiquidTube.Q_flow = -logDtHot*alphaHot*surfaceArea;
  heatFlowHotLiquidTube.Q_flow+heatFlowColdLiquidTube.Q_flow =0;
  //heatFlowColdLiquidTube.Q_flow = -logDtCold*alphaCold*surfaceArea;

   -alphaHot*(T_h_out - TwallCold) + der(TwallCold)*plateMetalMass*cpMetal/surfaceArea = -alphaCold*(TwallCold-T_c_in);
   -alphaHot*(T_h_in-TwallHot)+der(TwallHot)*plateMetalMass*cpMetal/surfaceArea  = -alphaCold*(TwallHot - T_c_out);

   QFlowFromHotLiqidOutToColdWall =  surfaceArea*alphaHot*(T_h_out - TwallCold);
   QFlowFromHotLiqidInToHotWall =  surfaceArea*alphaHot*(T_h_in-TwallHot);
   QFlowFromColdLiqidInToColdWall = - surfaceArea*alphaCold*(TwallCold- T_c_in);
   QFlowFromColdLiqidOutToHotWall =  -surfaceArea*alphaCold*(TwallHot - T_c_out);

  // connect(heatFlowHotLiquid,hotLiquidTube.heatPort);
   // connect(heatFlowColdLiquid,coldLiquidTube.heatPort);
   connect(alpha_in,alpha_internal);

    connect(hotLiquidIn, hotLiquidTube.liquidIn)
                                               annotation (Line(
        points={{70,70},{58,70},{58,30},{-54,30},{-54,19.6}},
        color={255,0,0},
        smooth=Smooth.None));
    connect(hotLiquidTube.liquidOut, hotLiquidOut)
                                                 annotation (Line(
        points={{-54,-19.6},{-54,-90},{-90,-90}},
        color={255,0,0},
        smooth=Smooth.None));
    connect(coldLiquidIn, coldLiquidTube.liquidIn)
                                                 annotation (Line(
        points={{-90,70},{-90,-40},{38,-40},{38,-19.6}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(coldLiquidTube.liquidOut, coldLiquidOut)
                                                   annotation (Line(
        points={{38,19.6},{38,24},{70,24},{70,-90}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(heatFlowColdLiquidTube.port, coldLiquidTube.heatPort) annotation (
        Line(
        points={{20,-2},{22,-2},{22,9.79717e-016},{30,9.79717e-016}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(heatFlowHotLiquidTube.port, hotLiquidTube.heatPort) annotation (Line(
        points={{-38,-1},{-36,-1},{-36,0},{-46,0}},
        color={191,0,0},
        smooth=Smooth.None));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                     graphics));
  end CcounterFlowHeatExchanger;

  model LiquidTube "ideal - no pressure loss"

    //import CCS.ThermalSeparation;
      outer ThermalSeparation.SystemTS systemTS;
  parameter Modelica.SIunits.Temperature T_ref=systemTS.T_ref "reference temperature"
                                                                          annotation(Dialog(tab="Advanced"));

  replaceable package MediumLiquid =
  ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O     constrainedby
      ThermalSeparation.Media.BaseMediumLiquid                                                annotation(choicesAllMatching);
    MediumLiquid.BaseProperties mediumLiquid(T0=T_ref, p=liquidIn.p, T=T_l, x=x_l,h=h_l);
      MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref, p=liquidIn.p, T=T_l_in, x=x_l_in,h=h_l_in);

    ThermalSeparation.Interfaces.LiquidPortIn liquidIn(redeclare package Medium =
          MediumLiquid) annotation (Placement(transformation(extent={{-10,88},{10,
              108}}, rotation=0), iconTransformation(extent={{-10,88},{10,108}})));

    ThermalSeparation.Interfaces.LiquidPortOut liquidOut(redeclare package Medium =
          MediumLiquid) annotation (Placement(transformation(extent={{-10,-108},{
              10,-88}}, rotation=0), iconTransformation(extent={{-10,-108},{10,-88}})));

  public
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort
      annotation (Placement(transformation(extent={{70,-24},{90,-4}}),
          iconTransformation(extent={{40,-10},{60,10}})));

  //parameter Boolean inertLiquid[nSL] = fill(false,nSL);

    final parameter Integer nSL = MediumLiquid.nSubstance;

  /*** Medium properties ***/
  Modelica.SIunits.Density rho_l;
  parameter Modelica.SIunits.Concentration c_l_start[nSL]={5,40000};
  Modelica.SIunits.Concentration c_l[nSL](start=c_l_start, each nominal=1e4);
  //Modelica.SIunits.Concentration c_l_in[nSL](start=c_l_start, each nominal=1e4);
  Modelica.SIunits.MoleFraction x_l[nSL];
    ThermalSeparation.Units.MolarEnthalpy h_l(stateSelect=StateSelect.prefer,start=1e6);
    ThermalSeparation.Units.MolarEnthalpy u_l;
    Modelica.SIunits.MolarMass MM_l;
  Modelica.SIunits.MoleFraction x_l_in[nSL];
    ThermalSeparation.Units.MolarEnthalpy h_l_in;

      Modelica.SIunits.Temperature T_l_in;
    Modelica.SIunits.Temperature T_l;
                       //(stateSelect=StateSelect.prefer);

    Modelica.SIunits.VolumeFlowRate Vdot_l(start=1e-3, nominal=1e3)=-liquidOut.Ndot*MM_l/rho_l;
    Modelica.SIunits.VolumeFlowRate Vdot_l_in=liquidIn.Ndot*mediumLiquidIn.MM/mediumLiquidIn.d;

    Modelica.SIunits.HeatFlowRate Qdot_wall;
    Modelica.SIunits.Pressure p_out;
    Modelica.SIunits.Pressure p_in;

    parameter Boolean ss_c2=false "true if c_l[2] is to be a state";
    Modelica.SIunits.Concentration dummy(stateSelect=StateSelect.prefer)=c_l[2] if
         ss_c2;
  //    SI.Concentration dummy2(stateSelect=StateSelect.prefer)=c_l[3] if ss_c2;

  public
    Real a = Qdot_wall/(2500*max(1e-6,Vdot_l_in)*1000);

    /*** geometry data ***/
      final parameter Modelica.SIunits.Area A=Modelica.Constants.pi/4*d_volume^2;
    final parameter Modelica.SIunits.Height H=length_out;
  parameter Real N_w= 20 "number of resistances in flow direction";
  //  parameter SI.Length s = 0.03 "distance between the center of two tubes";
    parameter Modelica.SIunits.Length d=0.01 "tube diameter";
   // final parameter Real a_pL = s/d;

    parameter Modelica.SIunits.Length height_in=-0.15;
    parameter Real zeta = 0.048;
    parameter Modelica.SIunits.Length length_out=0.15;
    parameter Modelica.SIunits.Diameter d_volume=0.025;
    // Real test;

  Real delta_h;
  Real delta_zeta;
  parameter Modelica.SIunits.Temperature T_start=300;
  parameter Boolean p_state=true "p is state";
  parameter Boolean useX_start = true;

  /*** Monitoring ***/
  Modelica.SIunits.Volume V_liq=A*H;
  Real sum_x = sum(x_l);

  equation
     //downstream
  //   liquidIn.p = liquidOut.p;
  liquidIn.x_outflow = inStream(liquidOut.x_outflow);
  liquidIn.h_outflow = inStream(liquidOut.h_outflow);

     liquidIn.p = p_in;
     liquidOut.p = p_out;

     liquidOut.h_outflow = h_l;
     inStream(liquidIn.x_outflow) = x_l_in;
     liquidOut.x_outflow = x_l;

  h_l_in=inStream(liquidIn.h_outflow);
    //h_l_in =mediumLiquidIn.h;

    //h_l =mediumLiquid.h;
    u_l =mediumLiquid.u;
     rho_l = mediumLiquid.d;
     MM_l = mediumLiquid.MM;

    /*** correlation between mole fraction x and concentration c ***/
    for i in 1:nSL loop
      x_l[i] = c_l[i] * MM_l /rho_l;
    end for;

      /***energy balance ***/
  A*H * der(sum(c_l[:])*h_l)  =liquidIn.Ndot*h_l_in + liquidOut.Ndot*h_l + Qdot_wall;

      /*** static mass balance ***/
      /*** if a dynamic mass balance is used, the concentration c_l needs initial values ***/
        liquidIn.Ndot + liquidOut.Ndot =0;
       x_l = x_l_in;
  //      sum(x_l)=1;
        //A*H*der(c_l) = Vdot_l_in*(c_l_in)- Vdot_l * (c_l);

    /*** heat port ***/
    Qdot_wall = heatPort.Q_flow;
    T_l=heatPort.T;

  //  p_out=p;

    /*** pressure loss ***/
    // Vdot_v=1;
  //p_in-p_out = 5* N_w *(eps_vap*rho_v+eps_liq*rho_l)/2*a_pL/(a_pL-1) * (Vdot_v+Vdot_l)/A;
  //p_in-p_out = height_in*Modelica.Constants.g_n*rho_l + zeta*length_out/(2*d*A^2) * (Vdot_l)^2 *  rho_l;
  p_in-p_out = 0;

   delta_h=height_in*Modelica.Constants.g_n*rho_l;
   delta_zeta=zeta*length_out/(2*d*A^2) * (Vdot_l)^2 *  rho_l;

  initial equation
    T_l=T_start;
    if useX_start then
      x_l[2]=1-1e-5;
   //  x_l[2]= 0.9446;
   //  x_l[3]=0.0157;
    end if;

    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                        graphics), Icon(coordinateSystem(preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics={Rectangle(
            extent={{-91,40},{91,-40}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,149,255},
            origin={1,0},
            rotation=90), Text(
            extent={{-67,28},{67,-28}},
            lineColor={0,0,0},
            textString="%name",
            origin={-5,6},
            rotation=90)}));
  end LiquidTube;
end HeatExchanger;
