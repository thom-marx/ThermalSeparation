within ThermalSeparation.PressureLoss.PlateColumn;
model Vdot_DryHydro
  "calculation of the volume flow rate using the dry pressure loss and the hydrostatic pressure loss"
 //"extends ..." wurde auskommentiert, damit man das Modell nicht auswählen kann, weil das nicht so dolle ist

 // extends BasicPressureLossPlate;
SI.Pressure deltaP_liq[n] "pressure loss due to the liquid on the tray";

equation
for j in 1:n loop

  deltaP_liq[j] = h[j] * eps_liq_2ph[j] * rho_l[j]* Modelica.Constants.g_n;
/*** equation from Stichlmair, Chem.Ing.Techn. 50 (1978) ***/
  Vdot[j] = sign(p[j]-p[j+1])*sqrt((abs(p[j]-p[j+1]) - deltaP_liq[j]) *2/rho_v[j]/ geometry.zeta * geometry.A_free^2);
end for;

  p_v_in = p[1] + Vdot_in*Vdot_in/geometry.A_free^2 * geometry.zeta * rho_v[1]/2 + h[1]*rho_l[n]*Modelica.Constants.g_n;
end Vdot_DryHydro;
