within ThermalSeparation.Components.Tanks;
model Tank2 "wie Tank, aber nur mit 2 Ports"

  Interfaces.LiquidPortIn portIn(nS=nSL) 
    annotation (Placement(transformation(extent={{-30,68},{-10,88}})));

replaceable package MediumLiquid = 
ThermalSeparation.Media.ModelicaMedia.Water.Absorption.CO2_H2O 
                                                     constrainedby
    ThermalSeparation.Media.BaseMediumLiquid                                                          annotation(choicesAllMatching);
   outer ThermalSeparation.SystemTS systemTS;
 parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));
 parameter Real f_psat = systemTS.f_psat "factor to scale p_sat water" annotation(Dialog(tab="Advanced"));
 parameter Boolean psat_Antoine= systemTS.psat_Antoine
    "true if p_sat of IF97 is used"                                              annotation(Dialog(tab="Advanced"));

 MediumLiquid.BaseProperties mediumLiquid(T0=T_ref, f_psat=f_psat, psat_Antoine = psat_Antoine,p=p, T=portOut.T, x=x_l, x_total=x_l);
    MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref, f_psat=f_psat, psat_Antoine = psat_Antoine,p=portIn.p, T=portIn.T, x=portIn.x, x_total=portIn.x);

    parameter Boolean inertLiquid[nSL] = fill(false,nSL);
  parameter Integer nS=2
    "number of species which are equal in vapour and liquid phase";
  parameter Integer nL=0
    "number of additional substances which are only in liquid phase";

  final parameter Integer nSL = nS + nL;

/*** Medium properties ***/
 SI.Density rho_l = mediumLiquid.d;
 SI.Density rho_l_in = mediumLiquidIn.d;
  SI.Concentration dummy(stateSelect=StateSelect.always)=c_l[1];
   ThermalSeparation.Units.MolarEnthalpy h_l(stateSelect=StateSelect.always)= mediumLiquid.h;
   ThermalSeparation.Units.MolarEnthalpy u_l = mediumLiquid.u;
   SI.MolarMass MM_l = mediumLiquid.MM;
   SI.MolarMass MM_l_in = mediumLiquidIn.MM;
        SI.Concentration c_l_in[nSL] = portIn.c;
   ThermalSeparation.Units.MolarEnthalpy h_l_in= mediumLiquidIn.h;

    SI.Concentration c_l[nSL](stateSelect=StateSelect.default);
    SI.MoleFraction x_l[nSL];
  SI.Pressure p(start=2e5);

  /*** geometry data ***/
    final parameter SI.Area A= Modelica.Constants.pi/4* d_volume^2;
  SI.Height level(stateSelect=StateSelect.always);

  parameter SI.Diameter d_volume = 0.025;
  parameter Real zeta=2;

  parameter SI.Pressure p_ambient = 1e5;

  Interfaces.LiquidPortOut portOut(nS=nSL) 
    annotation (Placement(transformation(extent={{-34,-78},{-14,-58}})));

  Modelica.Blocks.Interfaces.RealOutput y=level 
    annotation (Placement(transformation(extent={{74,58},{94,78}})));
/*** Initialization ***/
parameter Real x_l_start[nSL]={1e-5, 1-1e-5} annotation(Dialog(tab="Initialization"));
parameter SI.Temperature T_start=300 annotation(Dialog(tab="Initialization"));
parameter SI.Height level_start= 0.1 annotation(Dialog(tab="Initialization"));

/*** Monitoring ***/
SI.Volume V_liq = A*level;
Real sum_x = sum(x_l);
SI.MassFlowRate mdot_in = portIn.Vdot*1000;
SI.MassFlowRate mdot_out=portOut.Vdot*1000;

equation
  portOut.x = x_l;
  portOut.c = c_l;
  portOut.p = p;

  /*** correlation between mole fraction x and concentration c ***/
   for i in 1:nSL loop
     x_l[i] = c_l[i] * MM_l /rho_l;
   end for;

     /***energy balance ***/
   A* der(sum(c_l[:])*h_l*level)  =   portIn.Vdot*sum(c_l_in[:])*h_l_in  + portOut.Vdot*sum(c_l[:])*h_l;
    /*** mass balance ***/
    for i in 1:nSL loop
            A* der(c_l[i]*level) = portIn.Vdot*c_l_in[i] + portOut.Vdot*c_l[i];
       end for;

/*** fill level ***/
A*der(level*rho_l/MM_l) = portIn.Vdot*rho_l_in/MM_l_in  + portOut.Vdot*rho_l/MM_l;

p=level*9.8*rho_l+portIn.p - zeta* rho_l/2*(portIn.Vdot/A)^2;

initial equation
   portOut.T=T_start;
   x_l=x_l_start;
  level=level_start;

  annotation (Diagram(graphics));
end Tank2;
