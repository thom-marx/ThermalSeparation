within ThermalSeparation.Components.LiquidVolumes;
model Volume "Liquid volume, constant fill level"
 extends Icons.Icons.Tank;
  ThermalSeparation.Interfaces.LiquidPortIn
                          portIn(redeclare package Medium=MediumLiquid)
    annotation (Placement(transformation(extent={{100,-20},{124,4}}),
        iconTransformation(extent={{-20,76},{20,116}})));

replaceable package MediumLiquid =
ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O     constrainedby
    ThermalSeparation.Media.BaseMediumLiquid                                                          annotation(choicesAllMatching);
   outer ThermalSeparation.SystemTS systemTS;
 parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

 MediumLiquid.BaseProperties mediumLiquid(T0=T_ref,p=p, T=T_l_out, x=x_l);
    MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref, p=portIn.p, T=T_l_in, x=x_l_in);

  parameter Boolean inertLiquid[nSL] = fill(false,nSL);

  final parameter Integer nSL = MediumLiquid.nSubstance
    "number of substances in the liquid phase";

/*** Medium properties ***/
 SI.Density rho_l = mediumLiquid.d;
 // SI.Concentration dummy[nSL-1](stateSelect=StateSelect.always)={c_l[1],c_l[2]};
   ThermalSeparation.Units.MolarEnthalpy h_l(stateSelect=StateSelect.default)= mediumLiquid.h;
   ThermalSeparation.Units.MolarEnthalpy u_l(stateSelect=StateSelect.always) = mediumLiquid.u;
   SI.MolarMass MM_l = mediumLiquid.MM;
        SI.Concentration x_l_in[nSL] = inStream(portIn.x_outflow);
   ThermalSeparation.Units.MolarEnthalpy h_l_in= mediumLiquidIn.h;
   SI.Temperature T_l_out;
   SI.Temperature T_l_in;
   SI.Density rho_l_in = mediumLiquidIn.d;
   SI.MolarMass MM_l_in = mediumLiquidIn.MM;

    SI.Concentration c_l[nSL](stateSelect=StateSelect.default);
    SI.MoleFraction x_l[nSL](stateSelect=StateSelect.always);
  SI.Pressure p(start=2e5);

  /*** geometry data ***/
    final parameter SI.Area A= Modelica.Constants.pi/4* d_volume^2;
 parameter SI.Height level=10;

  parameter SI.Diameter d_volume = 0.025;
  parameter Real zeta=2;

  parameter SI.Pressure p_ambient = 1e5;

  ThermalSeparation.Interfaces.LiquidPortOut
                           portOut(redeclare package Medium=MediumLiquid)
    annotation (Placement(transformation(extent={{100,-100},{126,-74}}),
        iconTransformation(extent={{-20,-116},{20,-76}})));

  Modelica.Blocks.Interfaces.RealOutput y=level
    annotation (Placement(transformation(extent={{80,20},{100,40}}),
        iconTransformation(extent={{80,20},{100,40}})));
/*** Initialization ***/
parameter Real x_l_start[nSL]={1e-5, 1-1e-5}  annotation(Dialog(tab="Initialization"));
parameter SI.Temperature T_start=300 annotation(Dialog(tab="Initialization"));
//parameter SI.Height level_start= 1 annotation(Dialog(tab="Initialization"));

/*** Monitoring ***/
SI.Volume V_liq = A*level;
Real sum_x = sum(x_l);
SI.MassFlowRate mdot_in = portIn.Ndot*mediumLiquidIn.MM;
SI.MassFlowRate mdot_out=-portOut.Ndot*mediumLiquid.MM;
SI.Concentration sum_c = sum(c_l);

equation
  portIn.h_outflow = h_l;
  portIn.x_outflow = x_l;

  h_l_in = inStream(portIn.h_outflow);

  portOut.x_outflow = x_l;
  portOut.h_outflow = h_l;
  portOut.p = p;

  /*** correlation between mole fraction x and concentration c ***/
   for i in 1:nSL loop
     x_l[i] = c_l[i] * MM_l /rho_l;
   end for;

     /***energy balance ***/
    // portIn.T=portOut.T;
   A* der(sum(c_l[:])*u_l*level)  = portIn.Ndot*h_l_in  +portOut.Ndot*h_l;
    /*** substance mole balance balance ***/
    for i in 1:nSL loop
            A* der(c_l[i]*level) = portIn.Ndot*inStream(portIn.x_outflow[i]) + portOut.Ndot * portOut.x_outflow[i];
    end for;
   // portIn.x=portOut.x;

/*** total mole balance ***/
A*der(level*rho_l/MM_l)= portIn.Ndot + portOut.Ndot;
//sum(x_l)=1;
p=portIn.p;//level*9.8*rho_l+portIn.p - zeta* rho_l/2*(portIn.Vdot/A)^2;

initial equation
   T_l_out=T_start;
   for i in 1:nSL loop
   x_l[i]=x_l_start[i];
   end for;
  // level=level_start;

  annotation (Diagram(graphics), Icon(coordinateSystem(preserveAspectRatio=false,
                   extent={{-100,-100},{100,100}}),
                                      graphics));
end Volume;
