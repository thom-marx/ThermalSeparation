within ThermalSeparation.Components.HeatExchanger;
model LayerOfPlate "LayerOfPlate"
  import SI = Modelica.SIunits;

//  constant Real pi=Modelica.Constants.pi;
//  parameter Integer N(min=1)=2 "Number of nodes";
 // parameter Boolean linear_gradient=false "no logarithmic correction"
/*   parameter Boolean userdefinedmaterial=true 
    "define own fixed material properties" annotation (Dialog(group="Material"));  
  replaceable ThermoPower.Thermal.MaterialProperties.Metals.CarbonSteel_A106C[N] Material(
              each npol=3) 
   extends ThermoPower.Thermal.MaterialProperties.Interfaces.PartialMaterial 
    "pre-defined material properties"       annotation (choicesAllMatching = true, Dialog(enable=userdefinedmaterial==false, group="Material")); */
   parameter SiemensPower.Utilities.Structures.metal metal
    "Wall metal properties"                                                      annotation (Dialog(enable=userdefinedmaterial, group="Material"));
//  parameter Integer Nt(min=1)=1 "Number of parallel tubes";
  parameter SI.Length L=1 "length";
  parameter SI.Length w=0.08 "width";
  parameter SI.Length s=0.008 "Wall thickness";

//  parameter Boolean WallRes=true "Wall conduction resistance accounted for";
//  parameter Boolean WithAxialHeatTransfer=false
//    "With heat transfer in the wall parallel to the flow direction"
  //        annotation (Dialog(enable=WallRes));
  parameter String initOpt="steadyState" "Initialisation option" annotation (Dialog(group="Initialization"),
  choices(
    choice="noInit" "No initial equations",
    choice="steadyState" "Steady-state initialization",
    choice="fixedTemperature" "Fixed-temperatures initialization"));
  parameter SI.Temperature T0 "Temperature start values"      annotation (Dialog(group="Initialization"));

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_ext(T(start = T0))
    "Outer heat port" 
    annotation (Placement(transformation(extent={{-16,20},{16,48}}, rotation=0)));                               //(T(start = linspace(Tstart1,TstartN,N)))
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_int(T(start = T0))
    "Inner heat port" 
    annotation (Placement(transformation(extent={{-14,-48},{16,-20}}, rotation=
            0)));

 // SI.Area Am "Area of the metal tube cross-section";
  SI.Temperature T(start=T0) "Node temperatures";
 //   SI.Temperature T[N](start=Modelica.Media.Water.WaterIF97_ph.temperature_ph(0.5*(p0+p1),linspace(h0,h1,N)));
//  SI.Length rint;
 // SI.Length rext;
  SI.Mass mass;
  SI.HeatCapacity HeatCap "HeatCapacity of a Tube part";
//  SI.HeatFlowRate Q_flow_ax[N] "axial(parallel) heat transfer";
//  SI.Density rho[N];
//  SI.SpecificHeatCapacity cp[N] "Metal heat capacity per unit volume [J/kg.K]";
//  SI.ThermalConductivity lambda[N] "Thermal conductivity";

initial equation
  if initOpt == "noInit" then
 // nothing to do
  elseif initOpt == "steadyState" then
    der(T) = 0;
  elseif initOpt == "fixedTemperature" then // fixed temperatures at start
    T = T0;
  else
    assert(false, "Unsupported initialisation option");
  end if;

equation

 mass=(metal.rho*L*w*s);
 HeatCap=metal.cp*mass;

  HeatCap * der(T) = port_int.Q_flow + port_ext.Q_flow "Energy balance";
  port_int.Q_flow = L*w*metal.lambda/s/2*(port_int.T-T);
  port_ext.Q_flow = L*w*metal.lambda/s/2*(port_ext.T-T);

                                annotation (Dialog(enable=WallRes),
              Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={
        Rectangle(
          extent={{-80,20},{80,-20}},
          lineColor={0,0,0},
          fillColor={128,128,128},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-84,-22},{-32,-50}},
          lineColor={0,0,0},
          fillColor={128,128,128},
          fillPattern=FillPattern.Forward,
          textString="Int"),
        Text(
          extent={{-82,50},{-34,24}},
          lineColor={0,0,0},
          fillColor={128,128,128},
          fillPattern=FillPattern.Forward,
          textString="Ext"),
        Text(
          extent={{-100,-60},{100,-88}},
          lineColor={191,95,0},
          textString="%name")}),
                           Documentation(info="<HTML>
<p>This is the model of a plate layer of solid material.
<p>The heat capacity (which is lumped at the center of the tube thickness) is accounted for, as well as the thermal resistance due to the finite heat conduction coefficient. Longitudinal heat conduction is neglected.
<p><b>Modelling options</b></p>
<p>The following options are available to specify the valve flow coefficient in fully open conditions:
<ul>
<li><tt>WallRes = false</tt>: the thermal resistance of the tube wall is neglected.
<li><tt>WallRes = true</tt>: the thermal resistance of the tube wall is accounted for.
</ul>
</HTML>",
        revisions="<html>
<ul>
<li> December 2006, adapted to SiemensPower by Haiko Steuer
<li><i>30 May 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Initialisation support added.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>
"), DymolaStoredErrors);
end LayerOfPlate;
