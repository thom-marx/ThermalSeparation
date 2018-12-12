within ThermalSeparation.Components.HeatExchanger;
model counterFlowHeatExchangerH2O_Karin
     extends Icons.Icons.HeatExchanger;
  import SI = Modelica.SIunits;

replaceable package MediumLiquidCold = 
      Modelica.Media.Water.WaterIF97_ph 
  constrainedby Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching);
replaceable package MediumLiquidHot = 
      ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O 
  constrainedby ThermalSeparation.Media.BaseMediumLiquid annotation(choicesAllMatching);

// parameter
    // Geometry
  parameter SI.Mass plateMetalMass = 310 annotation(Dialog(group="Geometry"));
  parameter SI.SpecificHeatCapacity cpMetal = 500 annotation(Dialog(group="Geometry"));
  parameter SI.Volume volumeHotLiquid = 11e-3 annotation(Dialog(group="Geometry"));
  parameter SI.Volume volumeColdLiquid = 11e-3 annotation(Dialog(group="Geometry"));

  parameter SI.Area surfaceArea = 10.4 annotation(Dialog(group="Geometry"));

  parameter SI.CoefficientOfHeatTransfer alphaHot = 2000 annotation(Dialog(group="Geometry"));
  parameter SI.CoefficientOfHeatTransfer alphaCold = 2000 annotation(Dialog(group="Geometry"));

    // nominal operating point
  parameter SI.Pressure dpHotLiquid = 2e5 annotation(Dialog(group="nominal operating point"));
  parameter SI.Pressure dpColdLiquid = 2e5 annotation(Dialog(group="nominal operating point"));
  parameter SI.MassFlowRate mFlowHotLiquid_nom = 0.3 annotation(Dialog(group="nominal operating point"));
  parameter SI.MassFlowRate mFlowColdLiquid_nom = 0.3 annotation(Dialog(group="nominal operating point"));
  parameter SI.Temperature ThotLiquidOut_nom = ThotLiquid_start annotation(Dialog(group="nominal operating point"));
  parameter SI.Temperature TcoldLiquidOut_nom = TcoldLiquid_start annotation(Dialog(group="nominal operating point"));

  parameter SI.AbsolutePressure pLiquid_nom = 4e5 annotation(Dialog(group="nominal operating point"));

    // initial
  parameter SI.Temperature ThotLiquidIn_start = 104 + 273.15 annotation(Dialog(group="initial"));
  parameter SI.Temperature TcoldLiquidIn_start = 51 + 273.15 annotation(Dialog(group="initial"));
  parameter SI.Temperature ThotLiquid_start = 56 + 273.15 annotation(Dialog(group="initial"));
  parameter SI.Temperature TcoldLiquid_start = 98 + 273.15 annotation(Dialog(group="initial"));

    // Medium parameter
  parameter Integer nShot = 2
    "number of species which are equal in vapour and liquid phase";
  parameter Integer nScold = 2
    "number of species which are equal in vapour and liquid phase";

 // final parameter

  final parameter Integer nLhot = nSLhot - nShot
    "number of additional substances which are only in liquid phase";
  final parameter Integer nLcold = 0
    "number of additional substances which are only in liquid phase";
  final parameter Integer nSLhot = MediumLiquidHot.nSubstance;
  final parameter Integer nSLcold = nScold + nLcold;

  final parameter SI.Density rhoHotLiquid_nom = MediumLiquidHot.density_pT(pLiquid_nom, ThotLiquidOut_nom);
  final parameter SI.Density rhoColdLiquid_nom = MediumLiquidCold.density_pT(pLiquid_nom, TcoldLiquidOut_nom);

  final parameter SI.Length dummyDiameter = 0.0254 " diameter is abitrary"; // since dp and volume is adjusted using tube zeta resp. length
  final parameter SI.Length hotLiquidTubeLength = volumeHotLiquid/(0.25*Modelica.Constants.pi*dummyDiameter^2);
  final parameter SI.Length coldLiquidTubeLength = volumeColdLiquid/(0.25*Modelica.Constants.pi*dummyDiameter^2);
  final parameter SI.Area CrossSectionArea = Modelica.Constants.pi*dummyDiameter^2/4;
  final parameter Real zetaHotLiquid = (dpHotLiquid*CrossSectionArea^2*2*dummyDiameter)/((mFlowHotLiquid_nom/rhoHotLiquid_nom)^2*hotLiquidTubeLength*rhoHotLiquid_nom);
  final parameter Real zetaHotLiquid_haiko = Modelica.Constants.pi^2*dummyDiameter^5/(8*hotLiquidTubeLength)*(dpHotLiquid*rhoHotLiquid_nom/mFlowHotLiquid_nom^2);
  final parameter Real zetaColdLiquid = Modelica.Constants.pi^2*dummyDiameter^5/(8*coldLiquidTubeLength)*(dpColdLiquid*rhoColdLiquid_nom/mFlowColdLiquid_nom^2);

//variables
 SI.TemperatureDifference logDtHot;
  SI.TemperatureDifference logDtCold;
 SI.TemperatureDifference logDtHot_approx= ((hotLiquidOut_T - TwallCold) + (hotLiquidIn.T - TwallHot))/2;
 SI.TemperatureDifference logDtCold_approx=(coldLiquidOut.T - TwallHot)+(coldLiquidIn.T - TwallCold)/2;
// SI.TemperatureDifference logDtHot_ln = ((hotLiquidOut_T - TwallCold)-(hotLiquidIn.T - TwallHot))/Modelica.Math.log(min(1e-5,(hotLiquidOut_T - TwallCold)/(hotLiquidIn.T - TwallHot)));
// SI.TemperatureDifference logDtCold_ln= ((coldLiquidOut.T - TwallHot)-(coldLiquidIn.T - TwallCold))/Modelica.Math.log(min(1e-5,(coldLiquidOut.T - TwallHot)/(coldLiquidIn.T - TwallCold)));
 SI.HeatFlowRate Qdot_cold= -alphaCold*surfaceArea*(TwallCold-coldLiquidIn.T)+alphaHot*surfaceArea*(hotLiquidOut_T - TwallCold);
 SI.HeatFlowRate Qdot_cold1= -alphaCold*surfaceArea*(TwallCold-coldLiquidIn.T);
 SI.HeatFlowRate Qdot_cold2= alphaHot*surfaceArea*(hotLiquidOut_T - TwallCold);
 SI.HeatFlowRate Qdot_hot=-alphaCold*surfaceArea*(TwallHot - coldLiquidOut.T)+alphaHot*surfaceArea*(hotLiquidIn.T-TwallHot);
 SI.HeatFlowRate Qdot_hot1=-alphaCold*surfaceArea*(TwallHot - coldLiquidOut.T);
 SI.HeatFlowRate Qdot_hot2=alphaHot*surfaceArea*(hotLiquidIn.T-TwallHot);

//SI.TemperatureDifference logDt;
//SI.TemperatureDifference logDt_approx=((hotLiquidOut_T - coldLiquidIn.T) + (hotLiquidIn.T - coldLiquidOut.T))/2;
//SI.TemperatureDifference logDt_ln = ((hotLiquidOut_T - coldLiquidIn.T)-(hotLiquidIn.T - coldLiquidOut.T))/Modelica.Math.log(min(1e-5,(hotLiquidOut_T - coldLiquidIn.T)/(hotLiquidIn.T - coldLiquidOut.T)));

SI.Temperature TwallHot(start=(ThotLiquidIn_start+TcoldLiquid_start)/2);
SI.Temperature TwallCold(start=(TcoldLiquidIn_start+ThotLiquid_start)/2);

// components
  parameter SI.Concentration c_l_start_hot[nSLhot]={5,40000};

  SI.HeatFlowRate QFlowFromColdLiqidInToColdWall;
SI.HeatFlowRate QFlowFromColdLiqidOutToHotWall;
SI.HeatFlowRate QFlowFromHotLiqidInToHotWall;
SI.HeatFlowRate QFlowFromHotLiqidOutToColdWall;

  ThermalSeparation.Components.HeatExchanger.LiquidTube hotLiquidTube(
    redeclare replaceable package MediumLiquid = MediumLiquidHot,
    d = dummyDiameter,
    height_in = 0,
    zeta = zetaHotLiquid,
    length_out = hotLiquidTubeLength,
    d_volume = dummyDiameter,
    T_start = ThotLiquid_start,
    useX_start = false,
    c_l_start=c_l_start_hot) 
    annotation (Placement(transformation(extent={{-72,-20},{-40,20}})));
  Interfaces.LiquidPortIn hotLiquidIn(redeclare package Medium=MediumLiquidHot)
                                                                                 annotation (Placement(transformation(
          extent={{62,-72},{82,-52}}),  iconTransformation(extent={{62,-72},{82,
            -52}})));
  Interfaces.LiquidPortOut hotLiquidOut(redeclare package Medium = 
        MediumLiquidHot)                                                           annotation (Placement(transformation(
          extent={{-82,-72},{-62,-52}}),  iconTransformation(extent={{-82,-72},
            {-62,-52}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
    heatFlowColdLiquidTube 
    annotation (Placement(transformation(extent={{0,-10},{20,10}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
    heatFlowHotLiquidTube 
    annotation (Placement(transformation(extent={{-20,-10},{-38,8}})));

  Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = 
        MediumLiquidCold) 
    annotation (Placement(transformation(extent={{62,48},{82,68}}),
        iconTransformation(extent={{62,48},{82,68}})));
  Modelica.Fluid.Valves.ValveLinear valveLinear(
    redeclare package Medium = MediumLiquidCold,
    dp_nominal=dpColdLiquid,
    dp_start=dpColdLiquid,
    m_flow_start=mFlowColdLiquid_nom,
    m_flow_small=0.001,
    show_T=false,
    show_V_flow=false,
    m_flow_nominal=mFlowColdLiquid_nom,
    allowFlowReversal=true) 
    annotation (Placement(transformation(extent={{62,-20},{82,0}})));
  Modelica.Fluid.Vessels.ClosedVolume volume(
    use_portsData=false,
    use_HeatTransfer=true,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    massDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial,
    V=volumeColdLiquid,
    redeclare package Medium = MediumLiquidCold,
    p_start=pLiquid_nom,
    T_start=TcoldLiquid_start,
    m_flow_small=0.001,
    nPorts=2) 
    annotation (Placement(transformation(extent={{30,-10},{50,10}})));
  Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = 
        MediumLiquidCold) 
    annotation (Placement(transformation(extent={{-78,50},{-58,70}}),
        iconTransformation(extent={{-78,50},{-58,70}})));
  Modelica.Blocks.Sources.Constant const(k=1) 
    annotation (Placement(transformation(extent={{44,40},{56,52}})));
  Modelica.Fluid.Sensors.Temperature coldLiquidIn(redeclare package Medium = 
        MediumLiquidCold) 
    annotation (Placement(transformation(extent={{-64,66},{-44,86}})));
  Modelica.Fluid.Sensors.Temperature coldLiquidOut(redeclare package Medium = 
        MediumLiquidCold) 
    annotation (Placement(transformation(extent={{32,66},{52,86}})));
  inner Modelica.Fluid.System system 
    annotation (Placement(transformation(extent={{-4,-64},{16,-44}})));

Real hotLiquidOut_T;
initial equation
   TwallHot = (ThotLiquidIn_start+TcoldLiquid_start)/2;
   TwallCold = (TcoldLiquidIn_start+ThotLiquid_start)/2;

equation
  hotLiquidOut.T=hotLiquidOut_T;
//  if  (hotLiquidOut_T - TwallCold) < 0.5 or (hotLiquidIn.T - TwallHot) < 0.5 then
//    logDtHot = logDtHot_approx;
//  else
//    logDtHot = logDtHot_ln;
//  end if;
//
//  if (coldLiquidOut.T - TwallHot) > -0.5 or (coldLiquidIn.T - TwallCold) > -0.5 then
//      logDtCold =logDtCold_approx;
//  else
//    logDtCold = logDtCold_ln;
//  end if;

 if  (hotLiquidOut.T - TwallCold) < 0.5  or (hotLiquidIn.T - TwallHot) < 0.5 then
   logDtHot = ((hotLiquidOut.T - TwallCold) + (hotLiquidIn.T - TwallHot))/2;
 else
   logDtHot = 0;//((hotLiquidOut.T - TwallCold)-(hotLiquidIn.T - TwallHot))/Modelica.Math.log((hotLiquidOut.T - TwallCold)/(hotLiquidIn.T - TwallHot));
 end if;

 if (coldLiquidOut.T - TwallHot) > -0.5 or (coldLiquidIn.T - TwallCold) > -0.5 then
     logDtCold =(coldLiquidOut.T - TwallHot)+(coldLiquidIn.T - TwallCold)/2;
 else
   logDtCold = 0;//((coldLiquidOut.T - TwallHot)-(coldLiquidIn.T - TwallCold))/Modelica.Math.log((coldLiquidOut.T - TwallHot)/(coldLiquidIn.T - TwallCold));
 end if;

//   if (hotLiquidOut.T - coldLiquidIn.T) > -0.5 or (hotLiquidIn.T - coldLiquidIn.T) > -0.5 then
//      logDt =logDt_approx;
//  else
//    logDt = logDt_ln;
//  end if;

 heatFlowHotLiquidTube.Q_flow = -(QFlowFromHotLiqidOutToColdWall+QFlowFromHotLiqidInToHotWall);
heatFlowColdLiquidTube.Q_flow = -(QFlowFromColdLiqidInToColdWall+QFlowFromColdLiqidOutToHotWall);

//heatFlowHotLiquidTube.Q_flow = -logDt*alphaHot*surfaceArea;
//heatFlowColdLiquidTube.Q_flow = -heatFlowHotLiquidTube.Q_flow;

//TwallCold=TwallHot;
//der(TwallCold) = heatFlowHotLiquidTube.Q_flow+heatFlowColdLiquidTube.Q_flow;
 -alphaHot*(hotLiquidOut_T - TwallCold) + der(TwallCold)*plateMetalMass*cpMetal/surfaceArea = -alphaCold*(TwallCold-coldLiquidIn.T);
 -alphaHot*(hotLiquidIn.T-TwallHot)+der(TwallHot)*plateMetalMass*cpMetal/surfaceArea  = -alphaCold*(TwallHot - coldLiquidOut.T);

  QFlowFromHotLiqidOutToColdWall =  surfaceArea/2*alphaHot*(hotLiquidOut.T - TwallCold);
 QFlowFromHotLiqidInToHotWall =  surfaceArea/2*alphaHot*(hotLiquidIn.T-TwallHot);
 QFlowFromColdLiqidInToColdWall = - surfaceArea/2*alphaCold*(TwallCold-coldLiquidIn.T);
 QFlowFromColdLiqidOutToHotWall =  -surfaceArea/2*alphaCold*(TwallHot - coldLiquidOut.T);

  connect(hotLiquidIn, hotLiquidTube.liquidIn) 
                                             annotation (Line(
      points={{72,-62},{58,-62},{58,30},{-56,30},{-56,19.6}},
      color={255,0,0},
      smooth=Smooth.None));
  connect(hotLiquidTube.liquidOut, hotLiquidOut) 
                                               annotation (Line(
      points={{-56,-19.6},{-56,-62},{-72,-62}},
      color={255,0,0},
      smooth=Smooth.None));
  connect(hotLiquidOut, hotLiquidOut)   annotation (Line(
      points={{-72,-62},{-65,-62},{-65,-54},{-58,-54},{-58,-62},{-72,-62}},
      color={0,0,0},
      pattern=LinePattern.None,
      thickness=1,
      smooth=Smooth.None));
  connect(heatFlowHotLiquidTube.port, hotLiquidTube.heatPort) annotation (Line(
      points={{-38,-1},{-36,-1},{-36,0},{-48,0}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(const.y, valveLinear.opening) annotation (Line(
      points={{56.6,46},{72,46},{72,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(heatFlowColdLiquidTube.port, volume.heatPort) annotation (Line(
      points={{20,0},{30,0}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(valveLinear.port_b, port_b) annotation (Line(
      points={{82,-10},{84,-10},{84,58},{72,58}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(coldLiquidIn.port, port_a) 
                                    annotation (Line(
      points={{-54,66},{-78,66},{-78,60},{-68,60}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(coldLiquidOut.port, port_b) 
                                     annotation (Line(
      points={{42,66},{57,66},{57,58},{72,58}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(port_a, volume.ports[1]) annotation (Line(
      points={{-68,60},{-30,60},{-30,-10},{38,-10}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(volume.ports[2], valveLinear.port_a) annotation (Line(
      points={{42,-10},{62,-10}},
      color={0,127,255},
      smooth=Smooth.None));
//   annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
//             -100},{100,100}}),
//                       graphics), Icon(coordinateSystem(preserveAspectRatio=false,
//           extent={{-100,-100},{100,100}}), graphics={
//         Rectangle(
//           extent={{-30,80},{-50,-80}},
//           fillPattern=FillPattern.Solid,
//           fillColor={135,135,135},
//           pattern=LinePattern.None),
//         Rectangle(
//           extent={{10,80},{-10,-80}},
//           fillPattern=FillPattern.Solid,
//           fillColor={135,135,135},
//           pattern=LinePattern.None),
//         Rectangle(
//           extent={{48,80},{28,-80}},
//           fillPattern=FillPattern.Solid,
//           fillColor={135,135,135},
//           pattern=LinePattern.None),
//         Rectangle(
//           extent={{86,80},{66,-80}},
//           fillPattern=FillPattern.Solid,
//           fillColor={135,135,135},
//           pattern=LinePattern.None),
//         Rectangle(
//           extent={{-68,80},{-88,-80}},
//           fillPattern=FillPattern.Solid,
//           fillColor={135,135,135},
//           pattern=LinePattern.None),
//         Line(
//           points={{-100,68},{-100,86},{-60,86},{-60,-98},{20,-98},{20,86},{100,
//               86},{100,70}},
//           smooth=Smooth.None,
//           color={0,0,255},
//           thickness=0.5),
//         Line(
//           points={{100,-72},{100,-90},{60,-90},{60,94},{-20,94},{-20,-90},{-100,
//               -90},{-100,-74}},
//           smooth=Smooth.None,
//           color={255,0,0},
//           thickness=0.5)}));
  annotation (Icon(graphics));
end counterFlowHeatExchangerH2O_Karin;
