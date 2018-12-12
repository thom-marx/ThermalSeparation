within ThermalSeparation.Components.LiquidVolumes;
model Tank2Inlets "tank model with varying liquid level and two inlet ports"
 extends Icons.Icons.ExpansionTank;
  ThermalSeparation.Interfaces.LiquidPortIn
                          portIn(redeclare package Medium=MediumLiquid)
    annotation (Placement(transformation(extent={{100,-40},{120,-20}}),
        iconTransformation(extent={{-58,76},{-18,116}})));

replaceable package MediumLiquid =
ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O     constrainedby
    ThermalSeparation.Media.BaseMediumLiquid                                                          annotation(choicesAllMatching);
   outer ThermalSeparation.SystemTS systemTS;
   parameter SI.Pressure p_gas= 1e5
    "pressure exerted by the (inert) gas above the liquid";
 parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

 MediumLiquid.BaseProperties mediumLiquid(T0=T_ref, p=portOut.p, T=T_l, x=x_l);
    MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref,p=portIn.p, T=T_l_in, x=x_l_in);
       MediumLiquid.BaseProperties mediumLiquidInRecirc(T0=T_ref,p=portIn.p, T=T_l_in_recirc, x=x_l_in_recirc);

    parameter Boolean inertLiquid[nSL] = fill(false,nSL);
  parameter Integer nS=2
    "number of species which are equal in vapour and liquid phase";
 final parameter Integer nL=nSL-nS
    "number of additional substances which are only in liquid phase";

  final parameter Integer nSL = MediumLiquid.nSubstance;

/*** Medium properties ***/
  SI.Density rho_l = mediumLiquid.d;
  SI.Density rho_l_in = mediumLiquidIn.d;
  SI.Density rho_l_in_recirc = mediumLiquidInRecirc.d;
 SI.Concentration c_l[nSL](stateSelect=StateSelect.default);
  SI.Concentration dummy(stateSelect=StateSelect.always)=c_l[1];
 SI.MoleFraction x_l[nSL];
   ThermalSeparation.Units.MolarEnthalpy h_l(stateSelect=StateSelect.always)= mediumLiquid.h;
   ThermalSeparation.Units.MolarEnthalpy u_l = mediumLiquid.u;

 SI.MoleFraction x_l_in[nSL] = inStream(portIn.x_outflow);
 SI.MoleFraction x_l_in_recirc[nSL] = inStream(portInRecirc.x_outflow);

   SI.MolarMass MM_l = mediumLiquid.MM;
    SI.MolarMass MM_l_in = mediumLiquidIn.MM;
     SI.MolarMass MM_l_in_recirc = mediumLiquidInRecirc.MM;

   ThermalSeparation.Units.MolarEnthalpy h_l_in= mediumLiquidIn.h;
      ThermalSeparation.Units.MolarEnthalpy h_l_in_recirc= mediumLiquidInRecirc.h;

  SI.Pressure p(start=2e5);

  /*** geometry data ***/
    final parameter SI.Area A= Modelica.Constants.pi/4* d_volume^2;
  SI.Height level(stateSelect=StateSelect.always);

  parameter SI.Diameter d_volume = 0.025;
  parameter Real zeta=2;

  parameter SI.Pressure p_ambient = 1e5;

SI.MassFlowRate mdot_in = portIn.Ndot*MM_l_in;
SI.MassFlowRate mdot_in2=portInRecirc.Ndot*MM_l_in_recirc;
SI.MassFlowRate mdot_out=-portOut.Ndot*MM_l;
SI.VolumeFlowRate Vdot_out = -portOut.Ndot*MM_l/rho_l;
SI.VolumeFlowRate Vdot_in = portIn.Ndot*MM_l_in/rho_l_in;
SI.VolumeFlowRate Vdot_in_recirc = portInRecirc.Ndot*MM_l_in_recirc/rho_l_in_recirc;

SI.Temperature T_l_in;
SI.Temperature T_l_in_recirc;
SI.Temperature T_l;

  ThermalSeparation.Interfaces.LiquidPortOut
                           portOut(redeclare package Medium=MediumLiquid)
    annotation (Placement(transformation(extent={{100,-100},{120,-80}}),
        iconTransformation(extent={{-20,-116},{20,-76}})));
  ThermalSeparation.Interfaces.LiquidPortIn
                          portInRecirc(redeclare package Medium=MediumLiquid)
    annotation (Placement(transformation(extent={{100,-20},{120,0}}),
        iconTransformation(extent={{20,76},{60,116}})));

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

equation
  portIn.h_outflow = inStream(portOut.h_outflow);
  portIn.x_outflow = inStream(portOut.x_outflow);
    portInRecirc.h_outflow = inStream(portOut.h_outflow);
  portInRecirc.x_outflow = inStream(portOut.x_outflow);

  h_l_in = inStream(portIn.h_outflow);
   h_l_in_recirc = inStream(portInRecirc.h_outflow);
   portOut.h_outflow = h_l;
portOut.x_outflow=x_l;
portOut.p = p;
  /*** correlation between mole fraction x and concentration c ***/
   for i in 1:nSL loop
     x_l[i] = c_l[i] * MM_l /rho_l;
   end for;

     /***energy balance ***/
   A* der(sum(c_l[:])*h_l*level)  =   portIn.Ndot*h_l_in + portInRecirc.Ndot*h_l_in_recirc + portOut.Ndot*h_l;
    /*** amount of substance balance for each component ***/
    for i in 1:nSL loop
          // A* der(c_l[i]*level) = portIn.Vdot*c_l_in[i] + portInRecirc.Vdot*c_l_in_recirc[i] - Vdot_out*c_l[i];
          A*der(level*rho_l*x_l[i]/MM_l*MediumLiquid.MMX[i]) = portIn.Ndot * x_l_in[i]*MediumLiquid.MMX[i]  +portOut.Ndot*x_l[i]*MediumLiquid.MMX[i];
    end for;

/*** total amount of substance balance ***/
A*der(level*rho_l/MM_l) = portIn.Ndot + portInRecirc.Ndot +portOut.Ndot;

portOut.p=level*9.8*rho_l+portIn.p - zeta* rho_l/2*((Vdot_in + Vdot_in_recirc)/A)^2;
portInRecirc.p=portIn.p;
portIn.p =p_gas;

initial equation
   T_l=T_start;
   x_l=x_l_start;
  level=level_start;

  annotation (Diagram(graphics), Icon(coordinateSystem(preserveAspectRatio=false,
                   extent={{-100,-100},{100,100}}),
                                      graphics),
    Documentation(info="<html>
<p>This model is the same as <a href=\"Modelica://ThermalSeparation.Components.LiquidVolumes.Tank\">Tank</a>, but with two liquid inlet ports.</p>
</html>"));
end Tank2Inlets;
