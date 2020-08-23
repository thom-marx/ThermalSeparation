within ThermalSeparation.Components.Tanks;
model Tank

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

 MediumLiquid.BaseProperties mediumLiquid(T0=T_ref, f_psat=f_psat, psat_Antoine = psat_Antoine,p=portOut.p, T=portOut.T, x=x_l, x_total=x_l);
    MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref, f_psat=f_psat, psat_Antoine = psat_Antoine,p=portIn.p, T=portIn.T, x=portIn.x, x_total=portIn.x);
       MediumLiquid.BaseProperties mediumLiquidInRecirc(T0=T_ref, f_psat=f_psat, psat_Antoine = psat_Antoine,p=portIn.p, T=portInRecirc.T, x=portInRecirc.x, x_total=portInRecirc.x);

    parameter Boolean inertLiquid[nSL] = fill(false,nSL);
  parameter Integer nS=2
    "number of species which are equal in vapour and liquid phase";
  parameter Integer nL=0
    "number of additional substances which are only in liquid phase";

  final parameter Integer nSL = nS + nL;

/*** Medium properties ***/
  SI.Density rho_l = mediumLiquid.d;
  SI.Density rho_l_in = mediumLiquidIn.d;
  SI.Density rho_l_in_recirc = mediumLiquidInRecirc.d;
 SI.Concentration c_l[nSL](stateSelect=StateSelect.default);
  SI.Concentration dummy(stateSelect=StateSelect.prefer)=c_l[1];
 SI.MoleFraction x_l[nSL];
   ThermalSeparation.Units.MolarEnthalpy h_l(stateSelect=StateSelect.prefer)= mediumLiquid.h;
   ThermalSeparation.Units.MolarEnthalpy u_l = mediumLiquid.u;
   SI.Concentration c_l_in[nSL];
      SI.Concentration c_l_in_recirc[nSL];
// SI.MoleFraction x_l_in[nSL];
//   ThermalSeparation.Units.MolarEnthalpy h_l_in;
   SI.MolarMass MM_l = mediumLiquid.MM;
    SI.MolarMass MM_l_in = mediumLiquidIn.MM;
     SI.MolarMass MM_l_in_recirc = mediumLiquidInRecirc.MM;

   ThermalSeparation.Units.MolarEnthalpy h_l_in= mediumLiquidIn.h;
      ThermalSeparation.Units.MolarEnthalpy h_l_in_recirc= mediumLiquidInRecirc.h;

  SI.Pressure p(start=2e5);

  /*** geometry data ***/
    final parameter SI.Area A= Modelica.Constants.pi/4* d_volume^2;
  SI.Height level(stateSelect=StateSelect.prefer);

  parameter SI.Diameter d_volume = 0.025;
  parameter Real zeta=2;

  parameter SI.Pressure p_ambient = 1e5;

SI.MassFlowRate mdot_in = portIn.Vdot*1000;
SI.MassFlowRate mdot_in2=portInRecirc.Vdot*1000;
SI.MassFlowRate mdot_out=Vdot_out*1000;
SI.VolumeFlowRate Vdot_out;

  Interfaces.LiquidPortOut portOut(nS=nSL) 
    annotation (Placement(transformation(extent={{-34,-78},{-14,-58}})));
  Interfaces.LiquidPortIn portInRecirc(nS=nSL) 
    annotation (Placement(transformation(extent={{60,12},{80,32}})));

  Modelica.Blocks.Interfaces.RealOutput y=level 
    annotation (Placement(transformation(extent={{74,58},{94,78}})));
/*** Initialization ***/
parameter Real x_l_start[nSL]={1e-5, 1-1e-5} annotation(Dialog(tab="Initialization"));
parameter SI.Temperature T_start=300 annotation(Dialog(tab="Initialization"));
parameter SI.Height level_start= 0.1 annotation(Dialog(tab="Initialization"));

/*** Monitoring ***/
SI.Volume V_liq = A*level;
Real sum_x = sum(x_l);

equation
portOut.x=x_l;
portOut.c=c_l;
portOut.p = p;
portOut.Vdot = -Vdot_out;
  /*** correlation between mole fraction x and concentration c ***/
   for i in 1:nSL loop
     x_l[i] = c_l[i] * MM_l /rho_l;
     mediumLiquidIn.x[i] = c_l_in[i]*MM_l_in/rho_l_in;
     portInRecirc.x[i] = c_l_in_recirc[i]*MM_l_in_recirc/rho_l_in_recirc;
   end for;

     /***energy balance ***/
   A* der(sum(c_l[:])*h_l*level)  =   portIn.Vdot*sum(c_l_in[:])*h_l_in + portInRecirc.Vdot*sum(c_l_in_recirc[:])*h_l_in_recirc - Vdot_out*sum(c_l[:])*h_l;
    /*** amount of substance balance for each component ***/
    for i in 1:nSL loop
          // A* der(c_l[i]*level) = portIn.Vdot*c_l_in[i] + portInRecirc.Vdot*c_l_in_recirc[i] - Vdot_out*c_l[i];
           A*der(level*rho_l*x_l[i]/MM_l*MediumLiquid.MMX[i]) = portIn.Vdot * rho_l_in*portIn.x[i]/MM_l_in*MediumLiquid.MMX[i] - Vdot_out*rho_l*x_l[i]/MM_l*MediumLiquid.MMX[i];
    end for;

/*** total amount of substance balance ***/
A*der(level*rho_l/MM_l) = portIn.Vdot*rho_l_in/MM_l_in + portInRecirc.Vdot*rho_l_in_recirc/MM_l_in_recirc - Vdot_out*rho_l/MM_l;

portOut.p=level*9.8*rho_l+portIn.p - zeta* rho_l/2*((portIn.Vdot+portInRecirc.Vdot)/A)^2;
portInRecirc.p=1e5;

initial equation
   portOut.T=T_start;
   x_l=x_l_start;
  level=level_start;

  annotation (Diagram(graphics));
end Tank;
