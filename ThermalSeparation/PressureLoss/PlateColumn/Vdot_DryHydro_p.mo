within ThermalSeparation.PressureLoss.PlateColumn;
model Vdot_DryHydro_p
  "calculation of the volume flow rate using the dry pressure loss and the hydrostatic pressure loss2"

  extends BasicPressureLossPlate;
  SI.Pressure deltaP_liq[n] "pressure loss due to the liquid on the tray";
  SI.Density rho_v_av[n-1];

equation
// die max-Abfrage braucht man, für den Fall, daß man die Kolonne leer initialisiert und das Inertgas nicht mit abbildet: in dem Fall würde sonst ein Volumenstrom in die falsche Richtung entstehen

for j in 1:n loop
  deltaP_liq[j] = h[j] * eps_liq_2ph[j] * rho_l_calc[j]* Modelica.Constants.g_n;
end for;

for j in 1:n-1 loop
  rho_v_av[j]    = (rho_v_calc[j]+rho_v_calc[j+1])/2;
/*** equation from Stichlmair, Chem.Ing.Techn. 50 (1978) ***/
  Vdot[j] = sign(p[j]-p[j+1])*sqrt(max(1e-7,(abs(p[j]-p[j+1]) - (deltaP_liq[j]+deltaP_liq[j])/2) *2/rho_v_av[j]/ geometry.zeta * geometry.A_free^2));
end for;

/*** outlet: half of the n-th discreet element ***/
/*** equation from Stichlmair, Chem.Ing.Techn. 50 (1978) ***/
  Vdot[n] = sign(p[n]-p[n+1])*sqrt(max(1e-7,(abs(p[n]-p[n+1]) - deltaP_liq[n]/2) *2/rho_v_calc[n]/ (geometry.zeta/2) * geometry.A_free^2));

/*** entry: half of the first discreet element ***/
  p_v_in = p[1] + Vdot_in*Vdot_in/geometry.A_free^2 * geometry.zeta/2 * rho_v_calc[1]/2 + deltaP_liq[1]/2;

//   extends BasicPressureLossPlate;
//   SI.Pressure deltaP_liq[n] "pressure loss due to the liquid on the tray";
//   SI.Density rho_v_av[n-1];
//
// equation
// // die max-Abfrage braucht man, für den Fall, daß man die Kolonne leer initialisiert und das Inertgas nicht mit abbildet: in dem Fall würde sonst ein Volumenstrom in die falsche Richtung entstehen
//
// for j in 1:n loop
//   deltaP_liq[j] = h[j] * eps_liq_2ph[j] * rho_l_calc[j]* Modelica.Constants.g_n;
// end for;
//
// for j in 1:n-1 loop
//   rho_v_av[j]    = (rho_v_calc[j]+rho_v_calc[j+1])/2;
// /*** equation from Stichlmair, Chem.Ing.Techn. 50 (1978) ***/
//   Vdot[j] = sign(p[j]-p[j+1])*sqrt(max(1e-7,(abs(p[j]-p[j+1]) - (deltaP_liq[j]+deltaP_liq[j])/2) *2/rho_v_av[j]/ geometry.zeta * geometry.A_free^2));
// end for;
//
// /*** outlet: half of the n-th discreet element ***/
// /*** equation from Stichlmair, Chem.Ing.Techn. 50 (1978) ***/
//   Vdot[n] = sign(p[n]-p[n+1])*sqrt(max(1e-7,(abs(p[n]-p[n+1]) - deltaP_liq[n]/2) *2/rho_v_calc[n]/ (geometry.zeta/2) * geometry.A_free^2));
//
// /*** entry: half of the first discreet element ***/
//   p_v_in = p[1] + Vdot_in*Vdot_in/geometry.A_free^2 * geometry.zeta/2 * rho_v_calc[1]/2 + deltaP_liq[1]/2;
  annotation (Documentation(info="<html>
<p>The pressure loss correlation was taken from [1]. It uses a zeta-value which is specific to the plate type used and has to be supplied in the geometry record.</p>
<p><br/>[1] Stichlmair, Chem.Ing.Techn. 50 (1978)</p>
</html>"));
end Vdot_DryHydro_p;
