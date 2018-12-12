within ThermalSeparation.Components.Tanks;
model Volume "Liquid volume, constant fill level"

  Interfaces.LiquidPortIn portIn(nS=nSL) 
    annotation (Placement(transformation(extent={{-30,68},{-6,92}})));

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

  parameter Integer nSL = 2 "number of substances in the liquid phase";

/*** Medium properties ***/
 SI.Density rho_l = mediumLiquid.d;
 // SI.Concentration dummy[nSL-1](stateSelect=StateSelect.always)={c_l[1],c_l[2]};
   ThermalSeparation.Units.MolarEnthalpy h_l(stateSelect=StateSelect.default)= mediumLiquid.h;
   ThermalSeparation.Units.MolarEnthalpy u_l(stateSelect=StateSelect.always) = mediumLiquid.u;
   SI.MolarMass MM_l = mediumLiquid.MM;
        SI.Concentration c_l_in[nSL] = portIn.c;
   ThermalSeparation.Units.MolarEnthalpy h_l_in= mediumLiquidIn.h;
   SI.Density rho_l_in = mediumLiquidIn.d;
   SI.MolarMass MM_l_in = mediumLiquidIn.MM;

    SI.Concentration c_l[nSL](stateSelect=StateSelect.default);
    SI.MoleFraction x_l[nSL](stateSelect=StateSelect.always);
  SI.Pressure p(start=2e5);
  SI.VolumeFlowRate Vdot_l_out;

  /*** geometry data ***/
    final parameter SI.Area A= Modelica.Constants.pi/4* d_volume^2;
 parameter SI.Height level=10;

  parameter SI.Diameter d_volume = 0.025;
  parameter Real zeta=2;

  parameter SI.Pressure p_ambient = 1e5;

  Interfaces.LiquidPortOut portOut(nS=nSL) 
    annotation (Placement(transformation(extent={{-34,-78},{-8,-52}})));

  Modelica.Blocks.Interfaces.RealOutput y=level 
    annotation (Placement(transformation(extent={{74,58},{94,78}})));
/*** Initialization ***/
parameter Real x_l_start[nSL]={1e-5, 1-1e-5}  annotation(Dialog(tab="Initialization"));
parameter SI.Temperature T_start=300 annotation(Dialog(tab="Initialization"));
//parameter SI.Height level_start= 1 annotation(Dialog(tab="Initialization"));

/*** Monitoring ***/
SI.Volume V_liq = A*level;
Real sum_x = sum(x_l);
SI.MassFlowRate mdot_in = portIn.Vdot*mediumLiquidIn.d;
SI.MassFlowRate mdot_out=-portOut.Vdot*mediumLiquid.d;
SI.Concentration sum_c = sum(c_l);
SI.Concentration sum_c_in = sum(c_l_in);

equation
  portOut.x = x_l;
  portOut.c = c_l;
  portOut.p = p;
  portOut.Vdot = -Vdot_l_out;

  /*** correlation between mole fraction x and concentration c ***/
   for i in 1:nSL loop
     x_l[i] = c_l[i] * MM_l /rho_l;
   end for;

     /***energy balance ***/
    // portIn.T=portOut.T;
   A* der(sum(c_l[:])*u_l*level)  = portIn.Vdot*sum(c_l_in[:])*h_l_in  -Vdot_l_out *sum(c_l[:])*h_l;
    /*** substance mole balance balance ***/
    for i in 1:nSL loop
            A* der(c_l[i]*level) = portIn.Vdot*c_l_in[i] - Vdot_l_out * c_l[i];
    end for;
   // portIn.x=portOut.x;

/*** total mole balance ***/
//A*der(level)= portIn.Vdot - Vdot_l_out;
A*der(level*rho_l/MM_l)= portIn.Vdot*rho_l_in/MM_l_in - Vdot_l_out*rho_l/MM_l;
//sum(x_l)=1;
p=portIn.p;//level*9.8*rho_l+portIn.p - zeta* rho_l/2*(portIn.Vdot/A)^2;

initial equation
   portOut.T=T_start;
   for i in 1:nSL loop
   x_l[i]=x_l_start[i];
   end for;
  // level=level_start;

  annotation (Diagram(graphics));
end Volume;
