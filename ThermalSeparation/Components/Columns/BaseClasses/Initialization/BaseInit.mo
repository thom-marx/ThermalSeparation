within ThermalSeparation.Components.Columns.BaseClasses.Initialization;
partial model BaseInit "base class for initial equations"
replaceable package MediumLiquid =
    ThermalSeparation.Media.BaseMediumLiquid;
replaceable package MediumVapour = ThermalSeparation.Media.BaseMediumVapour;
parameter Integer nS;
parameter Integer mapping[nS,2];
parameter Integer n;
parameter Integer nSL;
parameter Boolean inertLiquid[nSL];
parameter Integer nSV;
parameter Boolean inertVapour[nSV];
 parameter Boolean considerStartUp;
  input MediumLiquid.ThermodynamicProperties propsLiq[n];
input SI.Pressure p_v[n];
input SI.Pressure p_v_start[n];
input Real Ndot_v_transfer[n,nSV];
input Real x_l[n,nSL];
input Real x_v[n,nSV];
input Real x_l_start[n,nSL];
input Real x_v_start[n,nSV];
input Real n_mol_L[n];
input Real n_mol_V[n];
input Real x_l_star[n,nSL];
input Real x_v_star[n,nSV];
input Real x_total_start[nSV];
input SI.Temperature T_v[n];
input SI.Temperature T_l[n];
input SI.Temperature T_v_start[n];
input SI.Temperature T_l_start[n];
input SI.Concentration c_l[n,nSL];
input Real Edot_l_transfer[n];
input Real rho_l[n];
input Real rho_v[n];

protected
Real v[n] = propsLiq.v;
equation

end BaseInit;
