within ThermalSeparation.Components.LiquidVolumes;
model Tank "tank model with varying liquid level"
 extends Icons.Icons.ExpansionTank;
  ThermalSeparation.Interfaces.LiquidPortIn
                          portIn(redeclare package Medium=MediumLiquid)
    annotation (Placement(transformation(extent={{100,0},{120,20}}),
        iconTransformation(extent={{-20,76},{20,116}})));

replaceable package MediumLiquid =
ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O     constrainedby
    ThermalSeparation.Media.BaseMediumLiquid "medium to be used"                                                         annotation(choicesAllMatching);
   outer ThermalSeparation.SystemTS systemTS;
      parameter SI.Pressure p_gas= 1e5
    "pressure exerted by the (inert) gas above the liquid";
 parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

 MediumLiquid.BaseProperties mediumLiquid(T0=T_ref,p=p, T=T_out, x=x_l,h=h_l);
    MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref, p=portIn.p, T=T_in, x=inStream(portIn.x_outflow),h=h_l_in);

    parameter Boolean inertLiquid[nSL] = fill(false,nSL)
    "true for each component which is inert";
  parameter Integer nS=2
    "number of species which are equal in vapour and liquid phase";
  final parameter Integer nL=nSL-nS
    "number of additional substances which are only in liquid phase";

  final parameter Integer nSL = MediumLiquid.nSubstance;

/*** Medium properties ***/
 SI.Density rho_l = mediumLiquid.d;
 SI.Density rho_l_in = mediumLiquidIn.d;
  SI.Concentration dummy(stateSelect=StateSelect.always)=c_l[1];
   ThermalSeparation.Units.MolarEnthalpy h_l(stateSelect=StateSelect.always);
   ThermalSeparation.Units.MolarEnthalpy u_l = mediumLiquid.u;
   SI.MolarMass MM_l = mediumLiquid.MM;
   SI.MolarMass MM_l_in = mediumLiquidIn.MM;
   SI.MoleFraction x_l_in[nSL] = inStream(portIn.x_outflow);
   ThermalSeparation.Units.MolarEnthalpy h_l_in= inStream(portIn.h_outflow);

    SI.Concentration c_l[nSL](each stateSelect=StateSelect.default);
  output SI.MoleFraction x_l[ nSL];
  SI.Pressure p(start=2e5);

  /*** geometry data ***/
    final parameter SI.Area A= Modelica.Constants.pi/4* d_volume^2;
  SI.Height level(stateSelect=StateSelect.always);

  parameter SI.Diameter d_volume = 0.025 "diameter of the tank";
  parameter Real zeta=2 "friction factor";

  parameter SI.Pressure p_ambient = 1e5 "pressure of ambience";
ResultsLiquidVolumes results(nSL = nSL, x_l=x_l, level=level, V_liq=V_liq, c_l = c_l, T=T_out) annotation(HideResult=false);
  ThermalSeparation.Interfaces.LiquidPortOut
                           portOut(redeclare package Medium=MediumLiquid)
    annotation (Placement(transformation(extent={{100,-100},{120,-80}}),
        iconTransformation(extent={{-20,-116},{20,-76}})));

  Modelica.Blocks.Interfaces.RealOutput y=level
    annotation (Placement(transformation(extent={{80,20},{100,40}}),
        iconTransformation(extent={{80,20},{100,40}})));
/*** Initialization ***/
parameter Real x_l_start[nSL]={1e-5, 1-1e-5} annotation(Dialog(tab="Initialization"));
parameter SI.Temperature T_start=300 annotation(Dialog(tab="Initialization"));
parameter SI.Height level_start= 0.1 annotation(Dialog(tab="Initialization"));

/*** Monitoring ***/
SI.Volume V_liq = A*level;
Real sum_x = sum(x_l);
SI.MassFlowRate mdot_in = portIn.Ndot*MM_l_in;
SI.MassFlowRate mdot_out=portOut.Ndot*MM_l;
SI.Temperature T_out;
SI.Temperature T_in;
SI.VolumeFlowRate Vdot_in = portIn.Ndot *MM_l_in/rho_l_in;

equation
  portIn.x_outflow = x_l;
  portIn.h_outflow = h_l;

  portOut.x_outflow = x_l;
  portOut.h_outflow = mediumLiquid.h;
  portOut.p = p;

//h_l_in=mediumLiquidIn.h;
  /*** correlation between mole fraction x and concentration c ***/
   for i in 1:nSL loop
     x_l[i] = c_l[i] * MM_l /rho_l;
   end for;

     /***energy balance ***/
   A* der(sum(c_l[:])*h_l*level)  =   portIn.Ndot*h_l_in  + portOut.Ndot*h_l;
    /*** mass balance ***/
     for i in 1:nSL loop
   A* der(c_l[i]*level) = portIn.Ndot*x_l_in[i]+ portOut.Ndot*x_l[i];
        end for;

/*** fill level ***/
A*der(level*rho_l/MM_l) = portIn.Ndot  + portOut.Ndot;

p=level*9.8*rho_l+portIn.p - zeta* rho_l/2*(Vdot_in/A)^2;

portIn.p =p_gas;

initial equation
   T_out=T_start;
   x_l=x_l_start;
  level=level_start;

  annotation (Diagram(graphics), Icon(coordinateSystem(preserveAspectRatio=false,
                   extent={{-100,-100},{100,100}}),
                                      graphics),
    Documentation(info="<html>
<p>This is a tank model for a liquid tank with a (inert) gas phase above the liquid level. This gas phase exerts a constant pressure on the liquid surface. The gas phase is not modelled explicitly and the pressure exerted by the gas phase is set by the user using the parameter p_gas. This pressure is then also propaged to the liquid inlet port. </p>
</html>"));
end Tank;
