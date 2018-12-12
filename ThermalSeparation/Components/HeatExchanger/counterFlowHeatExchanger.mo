within ThermalSeparation.Components.HeatExchanger;
model counterFlowHeatExchanger
     extends Icons.Icons.HeatExchanger;
  import SI = Modelica.SIunits;

replaceable package MediumLiquidCold = 
      ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O 
  constrainedby ThermalSeparation.Media.BaseMediumLiquid annotation(choicesAllMatching);
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
   SI.CoefficientOfHeatTransfer alphaHot_nom = 2000 annotation(Dialog(group="Geometry"));
   SI.CoefficientOfHeatTransfer alphaCold_nom = 2000 annotation(Dialog(group="Geometry"));

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
  final parameter Integer nShot = 2
    "number of species which are equal in vapour and liquid phase";
  parameter Integer nScold = 2
    "number of species which are equal in vapour and liquid phase";

 // final parameter

  final parameter Integer nLhot = nSLhot - nShot
    "number of additional substances which are only in liquid phase";
  final parameter Integer nLcold = nSLcold - nScold
    "number of additional substances which are only in liquid phase";
  final parameter Integer nSLhot = MediumLiquidHot.nSubstance;
  final parameter Integer nSLcold = MediumLiquidCold.nSubstance;

  final parameter SI.Density rhoHotLiquid_nom = MediumLiquidHot.density_pT(pLiquid_nom, ThotLiquidOut_nom);
  final parameter SI.Density rhoColdLiquid_nom = MediumLiquidCold.density_pT(pLiquid_nom, TcoldLiquidOut_nom);

  final parameter SI.Length dummyDiameter = 0.0254 " diameter is abitrary"; // since dp and volume is adjusted using tube zeta resp. length
  final parameter SI.Length hotLiquidTubeLength = volumeHotLiquid/(0.25*Modelica.Constants.pi*dummyDiameter^2);
  final parameter SI.Length coldLiquidTubeLength = volumeColdLiquid/(0.25*Modelica.Constants.pi*dummyDiameter^2);
  final parameter SI.Area CrossSectionArea = Modelica.Constants.pi*dummyDiameter^2/4;
  final parameter Real zetaHotLiquid = (dpHotLiquid*CrossSectionArea^2*2*dummyDiameter)/((mFlowHotLiquid_nom/rhoHotLiquid_nom)^2*hotLiquidTubeLength*rhoHotLiquid_nom);
  final parameter Real zetaHotLiquid_haiko = Modelica.Constants.pi^2*dummyDiameter^5/(8*hotLiquidTubeLength)*(dpHotLiquid*rhoHotLiquid_nom/mFlowHotLiquid_nom^2);
  final parameter Real zetaColdLiquid = Modelica.Constants.pi^2*dummyDiameter^5/(8*coldLiquidTubeLength)*(dpColdLiquid*rhoColdLiquid_nom/mFlowColdLiquid_nom^2);
//parameter SI.Area areaHT = 1 "heat transfer area";

  parameter SI.Concentration c_l_start_hot[nSLhot]={5,40000};
  parameter SI.Concentration c_l_start_cold[nSLcold]={5,40000};

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
parameter Boolean ss_hot_c2=true
    "true if c_l[2] of the hot liquid is to be a state";
parameter Boolean ss_cold_c2=true
    "true if c_l[2] of the cold liquid is to be a state";

// components

  ThermalSeparation.Components.HeatExchanger.LiquidTube hotLiquidTube(
    redeclare replaceable package MediumLiquid = MediumLiquidHot,
    d = dummyDiameter,
    height_in = 0,
    zeta = zetaHotLiquid,
    length_out = hotLiquidTubeLength,
    d_volume = dummyDiameter,
    T_start = ThotLiquid_start,
    useX_start = false,
    ss_c2=ss_hot_c2,
    c_l_start=c_l_start_hot) 
    annotation (Placement(transformation(extent={{-70,-20},{-38,20}})));
  ThermalSeparation.Components.HeatExchanger.LiquidTube coldLiquidTube(
    redeclare replaceable package MediumLiquid = MediumLiquidCold,
    d = dummyDiameter,
    height_in = 0,
    zeta = zetaColdLiquid,
    length_out = coldLiquidTubeLength,
    d_volume = dummyDiameter,
    T_start = TcoldLiquid_start,
    useX_start = false,
    ss_c2=ss_cold_c2,
    c_l_start=c_l_start_cold)  annotation (Placement(transformation(
        extent={{-16,-20},{16,20}},
        rotation=180,
        origin={38,0})));
  Interfaces.LiquidPortIn hotLiquidIn(redeclare package Medium=MediumLiquidHot)
                                                                                 annotation (Placement(transformation(
          extent={{60,-72},{80,-52}}),  iconTransformation(extent={{60,-72},{80,
            -52}})));
  Interfaces.LiquidPortOut hotLiquidOut(redeclare package Medium = 
        MediumLiquidHot)                                                           annotation (Placement(transformation(
          extent={{-82,-72},{-62,-52}}),  iconTransformation(extent={{-82,-72},
            {-62,-52}})));
  Interfaces.LiquidPortIn coldLiquidIn(redeclare package Medium = 
        MediumLiquidCold)                                                         annotation (Placement(transformation(
          extent={{-82,54},{-62,74}}),  iconTransformation(extent={{-82,54},{
            -62,74}})));
  Interfaces.LiquidPortOut coldLiquidOut( redeclare package Medium = 
        MediumLiquidCold)                                                             annotation (Placement(transformation(
          extent={{62,52},{82,72}}),  iconTransformation(extent={{62,52},{82,72}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
    heatFlowColdLiquidTube 
    annotation (Placement(transformation(extent={{0,-10},{20,10}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
    heatFlowHotLiquidTube 
    annotation (Placement(transformation(extent={{-20,-10},{-38,8}})));

initial equation
   TwallHot = (ThotLiquidIn_start+TcoldLiquid_start)/2;
   TwallCold = (TcoldLiquidIn_start+ThotLiquid_start)/2;

equation

alphaCold  = alphaCold_nom*(rhoColdLiquid_nom*coldLiquidIn.Vdot)/mFlowColdLiquid_nom;
alphaHot = alphaHot_nom*(rhoHotLiquid_nom*hotLiquidIn.Vdot)/mFlowHotLiquid_nom;
 if  (hotLiquidOut.T - TwallCold) < 0.1  or (hotLiquidIn.T - TwallHot) < 0.1 then
   logDtHot = ((hotLiquidOut.T - TwallCold) + (hotLiquidIn.T - TwallHot))/2;
 else
   logDtHot = ((hotLiquidOut.T - TwallCold)-(hotLiquidIn.T - TwallHot))/Modelica.Math.log((hotLiquidOut.T - TwallCold)/(hotLiquidIn.T - TwallHot));
 end if;

 if (coldLiquidOut.T - TwallHot) > -0.1 or (coldLiquidIn.T - TwallCold) > -0.1 then
     logDtCold =(coldLiquidOut.T - TwallHot)+(coldLiquidIn.T - TwallCold)/2;
 else
   logDtCold = ((coldLiquidOut.T - TwallHot)-(coldLiquidIn.T - TwallCold))/Modelica.Math.log((coldLiquidOut.T - TwallHot)/(coldLiquidIn.T - TwallCold));
 end if;

heatFlowHotLiquidTube.Q_flow = -logDtHot*alphaHot*surfaceArea;
heatFlowColdLiquidTube.Q_flow = -logDtCold*alphaCold*surfaceArea;

 -alphaHot*(hotLiquidOut.T - TwallCold) + der(TwallCold)*plateMetalMass*cpMetal/surfaceArea = -alphaCold*(TwallCold-coldLiquidIn.T);
 -alphaHot*(hotLiquidIn.T-TwallHot)+der(TwallHot)*plateMetalMass*cpMetal/surfaceArea  = -alphaCold*(TwallHot - coldLiquidOut.T);

 QFlowFromHotLiqidOutToColdWall =  surfaceArea*alphaHot*(hotLiquidOut.T - TwallCold);
 QFlowFromHotLiqidInToHotWall =  surfaceArea*alphaHot*(hotLiquidIn.T-TwallHot);
 QFlowFromColdLiqidInToColdWall = - surfaceArea*alphaCold*(TwallCold-coldLiquidIn.T);
 QFlowFromColdLiqidOutToHotWall =  -surfaceArea*alphaCold*(TwallHot - coldLiquidOut.T);

// connect(heatFlowHotLiquid,hotLiquidTube.heatPort);
 // connect(heatFlowColdLiquid,coldLiquidTube.heatPort);

  connect(hotLiquidIn, hotLiquidTube.liquidIn) 
                                             annotation (Line(
      points={{70,-62},{58,-62},{58,30},{-54,30},{-54,19.6}},
      color={255,0,0},
      smooth=Smooth.None));
  connect(hotLiquidTube.liquidOut, hotLiquidOut) 
                                               annotation (Line(
      points={{-54,-19.6},{-54,-62},{-72,-62}},
      color={255,0,0},
      smooth=Smooth.None));
  connect(coldLiquidIn, coldLiquidTube.liquidIn) 
                                               annotation (Line(
      points={{-72,64},{-70,64},{-70,-52},{38,-52},{38,-19.6}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(coldLiquidTube.liquidOut, coldLiquidOut) 
                                                 annotation (Line(
      points={{38,19.6},{38,62},{72,62}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(coldLiquidOut, coldLiquidOut)   annotation (Line(
      points={{72,62},{65.5,62},{65.5,54},{59,54},{59,62},{72,62}},
      color={0,0,0},
      pattern=LinePattern.None,
      thickness=1,
      smooth=Smooth.None));
  connect(hotLiquidOut, hotLiquidOut)   annotation (Line(
      points={{-72,-62},{-65,-62},{-65,-54},{-58,-54},{-58,-62},{-72,-62}},
      color={0,0,0},
      pattern=LinePattern.None,
      thickness=1,
      smooth=Smooth.None));
  connect(heatFlowColdLiquidTube.port, coldLiquidTube.heatPort) annotation (
      Line(
      points={{20,0},{22,0},{22,9.79717e-016},{30,9.79717e-016}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(heatFlowHotLiquidTube.port, hotLiquidTube.heatPort) annotation (Line(
      points={{-38,-1},{-36,-1},{-36,0},{-46,0}},
      color={191,0,0},
      smooth=Smooth.None));
//   annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
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
end counterFlowHeatExchanger;
