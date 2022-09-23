within ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy;
partial model BaseHomotopy

 parameter Integer nS=2 "number of substances crossing phase boundary" annotation(Dialog(tab="Propagated by column model",enable=false));
 parameter Integer nSL=2 "number of liquid substances" annotation(Dialog(tab="Propagated by column model",enable=false));
 parameter Integer nSV=2 "number of vapour substances"
                                                      annotation(Dialog(tab="Propagated by column model",enable=false));
 parameter Integer n=1 "number of stages" annotation(Dialog(tab="Propagated by column model",enable=false));
 Boolean useHomotopy=false;

 Modelica.Units.SI.MolarFlowRate Ndot_l_inter[n,nSL]=fill(0,n,nSL)
    "nominal molar flow across liquid phase boundary";
 Modelica.Units.SI.MolarFlowRate Ndot_v_inter[n,nSV]=fill(0,n,nSV)
    "nominal molar flow across vapour phase boundary";
 Modelica.Units.SI.HeatFlowRate Edot_l_inter[n]=fill(0,n)
    "nominal heat flow rate across liquid phase boundary";
 Modelica.Units.SI.HeatFlowRate Edot_v_inter[n]=fill(0,n)
    "nominal heat flow rate across vapour phase boundary";
 Modelica.Units.SI.Pressure dp=0.1e5 "nominal pressure loss";
 Real K[nS]=fill(1,nS) "nominal equilibrium constant";
 Modelica.Units.SI.Density rho_liq=1000 "nominal liquid density";
 Modelica.Units.SI.Density rho_vap=1 "nominal vapour density";
 ThermalSeparation.Units.MolarEnthalpy h_liq=2e3 "nominal liquid emthalpy";
 ThermalSeparation.Units.MolarEnthalpy h_vap=1.5e4 "nominal vapour enthalpy";

end BaseHomotopy;
