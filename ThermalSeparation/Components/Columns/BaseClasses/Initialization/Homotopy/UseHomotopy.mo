within ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy;
model UseHomotopy

extends
    ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.BaseHomotopy(
    useHomotopy=true,
    Ndot_l_inter=Ndot_l_inter_nom,
    Ndot_v_inter=Ndot_v_inter_nom,
    Edot_l_inter=Edot_l_inter_nom,
    Edot_v_inter=Edot_v_inter_nom,
    dp=dp_nom,
    K=K_nom);

 parameter Boolean whichHomotopy[6]=fill(false,6)
    "set true to use Homotopy on 1: Ndot_inter, 2:Edot_inter, 3: dp, 4: K, 5:rho, 6:h";
Boolean bool_Ndot_inter=if not whichHomotopy[1] then false else true
    "true if homotopy is applied on molar flow rate across phase boundary";
Boolean bool_Edot_inter=if not whichHomotopy[2] then false else true
    "true if homotopy is applied on heat flow rate across phase boundary";
Boolean bool_dp=if not whichHomotopy[3] then false else true
    "true if homotopy is applied on pressure loss";
Boolean bool_K=if not whichHomotopy[4] then false else true
    "true is homotopy is applied on equilibrium constant";
Boolean bool_rho=if not whichHomotopy[5] then false else true "true is homotopy is applied on liquid and vapour density";
Boolean bool_h=if not whichHomotopy[6] then false else true "true is homotopy is applied on liquid and vapour enthalpy";

parameter Modelica.SIunits.MolarFlowRate Ndot_l_inter_nom[n,nSL]=fill(0,n,nSL)
    "nominal molar flow across liquid phase boundary"                                                                        annotation(Dialog(group="Nominal values - used if respective whichHomotopy is set true"));
parameter Modelica.SIunits.MolarFlowRate Ndot_v_inter_nom[n,nSV]=fill(0,n,nSV)
    "nominal molar flow across vapour phase boundary"                                                                         annotation(Dialog(group="Nominal values - used if respective whichHomotopy is set true"));
parameter Modelica.SIunits.HeatFlowRate Edot_l_inter_nom[n]=fill(0,n)
    "nominal heat flow rate across liquid phase boundary"                                                                annotation(Dialog(group="Nominal values - used if respective whichHomotopy is set true"));
parameter Modelica.SIunits.HeatFlowRate Edot_v_inter_nom[n]=fill(0,n)
    "nominal heat flow rate across vapour phase boundary"                                                                annotation(Dialog(group="Nominal values - used if respective whichHomotopy is set true"));
parameter Modelica.SIunits.Pressure dp_nom=0.1e5 "nominal pressure loss" annotation(Dialog(group="Nominal values - used if respective whichHomotopy is set true"));
parameter Real K_nom[nS]=fill(1,nS) "nominal equilibrium constant" annotation(Dialog(group="Nominal values - used if respective whichHomotopy is set true"));
parameter Modelica.SIunits.Density rho_liq_nom=1000 "nominal liquid density"              annotation(Dialog(group="Nominal values - used if respective whichHomotopy is set true"));
parameter Modelica.SIunits.Density rho_vap_nom=1 "nominal vapour density"          annotation(Dialog(group="Nominal values - used if respective whichHomotopy is set true"));
parameter ThermalSeparation.Units.MolarEnthalpy h_liq_nom=2e3 "nominal liquid emthalpy"              annotation(Dialog(group="Nominal values - used if respective whichHomotopy is set true"));
parameter ThermalSeparation.Units.MolarEnthalpy h_vap_nom=1.5e4 "nominal vapour enthalpy"          annotation(Dialog(group="Nominal values - used if respective whichHomotopy is set true"));

end UseHomotopy;
