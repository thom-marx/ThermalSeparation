within ThermalSeparation.PressureLoss;
package TrayColumn
  partial model BasicPressureLossPlate
    "basic pressure loss model for a plate column"
    extends BasicPressureLoss;
      replaceable record Geometry =
        ThermalSeparation.Geometry.PlateColumn.Geometry                            constrainedby
      ThermalSeparation.Geometry.PlateColumn.Geometry annotation(Dialog(enable=false));
      Geometry geometry(n=n);

    input SI.Height h[n](start=0.9*ones(n))
      "height of the 2ph region on the tray";
    input Real eps_liq_2ph[n]
      "liquid fraction in the two-phase area on the tray";
    input Boolean startUp[n];
  equation

  end BasicPressureLossPlate;

  model DryHydrostatic
    "calculation of the volume flow rate using the dry pressure loss and the hydrostatic pressure loss2"

     extends BasicPressureLossPlate;
     SI.Pressure deltaP_liq[n] "pressure loss due to the liquid on the tray";
     SI.Density rho_v_av[n-1];

  equation
   for j in 1:n loop
     deltaP_liq[j] = h[j] * eps_liq_2ph[j] * rho_l_calc[j]* Modelica.Constants.g_n;
   end for;

   for j in 1:n-1 loop
     rho_v_av[j]    = (rho_v_calc[j]+rho_v_calc[j+1])/2;
   /*** equation from Stichlmair, Chem.Ing.Techn. 50 (1978), Signatur: Sti78 ***/
   end for;

  if homotopyMethod.bool_dp and homotopyMethod.useHomotopy then
   for j in 1:n-1 loop
  sign(p[j]-p[j+1])*max(1e-7,(abs(p[j]-p[j+1])))=homotopy(actual=Vdot[j]^2*sqrt(rho_v_av[j])* geometry.zeta / geometry.A_free^2+(deltaP_liq[j]+deltaP_liq[j])/2,simplified=homotopyMethod.dp);
   end for;

   /*** outlet: half of the n-th discreet element ***/
   /*** equation from Stichlmair, Chem.Ing.Techn. 50 (1978), Signatur: Sti78 ***/
  sign(p[n]-p[n+1])*max(1e-7,(abs(p[n]-p[n+1])))=homotopy(actual=Vdot[n]^2*sqrt(rho_v_calc[n])* (geometry.zeta/2) / geometry.A_free^2+deltaP_liq[n]/2,simplified=homotopyMethod.dp/2);

   /*** entry: half of the first discreet element ***/
     p_v_in-p[1] =  homotopy(actual=Vdot_in*Vdot_in/geometry.A_free^2 * geometry.zeta/2 * sqrt(rho_v_calc[1]) + deltaP_liq[1]/2,simplified=homotopyMethod.dp/2);

  else
   for j in 1:n-1 loop
  sign(p[j]-p[j+1])*max(1e-7,(abs(p[j]-p[j+1])))=Vdot[j]^2*sqrt(rho_v_av[j])* geometry.zeta / geometry.A_free^2+(deltaP_liq[j]+deltaP_liq[j])/2;

   end for;

   /*** outlet: half of the n-th discreet element ***/
   /*** equation from Stichlmair, Chem.Ing.Techn. 50 (1978), Signatur: Sti78 ***/
  sign(p[n]-p[n+1])*max(1e-7,(abs(p[n]-p[n+1])))=Vdot[n]^2*sqrt(rho_v_calc[n])* (geometry.zeta/2) / geometry.A_free^2+deltaP_liq[n]/2;

   /*** entry: half of the first discreet element ***/
     p_v_in-p[1] =  Vdot_in*Vdot_in/geometry.A_free^2 * geometry.zeta/2 * sqrt(rho_v_calc[1]) + deltaP_liq[1]/2;
  end if;
    annotation (Documentation(info="<html>
<p>The pressure loss correlation was taken from [1]. It uses a zeta-value which is specific to the plate type used and has to be supplied in the geometry record.</p>
<p><br/>[1] Stichlmair, Chem.Ing.Techn. 50 (1978)</p>
</html>"));
  end DryHydrostatic;

  model NominalLinear "linear pressure drop with nominal values"

    extends BasicPressureLossPlate;
    parameter SI.Pressure deltaP_nom=0.05e5 "nominal pressure drop";
    parameter SI.VolumeFlowRate Vdot_nom=0.4 "nominal volume flow rate";
    final parameter Real K=(deltaP_nom/Vdot_nom)/n;

  equation
    // die max-Abfrage braucht man, für den Fall, daß man die Kolonne leer initialisiert und das Inertgas nicht mit abbildet: in dem Fall würde sonst ein Volumenstrom in die falsche Richtung entstehen

  if homotopyMethod.bool_dp and homotopyMethod.useHomotopy then
    for j in 1:n - 1 loop
       max(0,sign(p[j] - p[j + 1])*max(0, (abs(p[j] - p[j + 1]))))=homotopy(actual=K*Vdot[j],simplified=homotopyMethod.dp);
    end for;

    /*** outlet: half of the n-th discreet element ***/
    max(0,sign(p[n] - p[n + 1])*max(0, (abs(p[n] - p[n + 1]))))=homotopy(actual=(K/2)*Vdot[n],simplified=homotopyMethod.dp/2);

    /*** entry: half of the first discreet element ***/
    p_v_in-p[1] = homotopy(actual=Vdot_in*K/2,simplified=homotopyMethod.dp/2);

  else
    for j in 1:n - 1 loop
       max(0,sign(p[j] - p[j + 1])*max(0, (abs(p[j] - p[j + 1]))))=K*Vdot[j];
    end for;

    /*** outlet: half of the n-th discreet element ***/
    max(0,sign(p[n] - p[n + 1])*max(0, (abs(p[n] - p[n + 1]))))=(K/2)*Vdot[n];

    /*** entry: half of the first discreet element ***/
    p_v_in-p[1] = Vdot_in*K/2;
  end if;

    annotation (Documentation(info="<html>
<p>Linear pressure drop, using a nominal value for delta_P and Vdot </p>
</html>"));
  end NominalLinear;

  model NominalLinear2 "linear pressure drop with nominal values2"

    extends BasicPressureLossPlate;
    parameter SI.Pressure deltaP_nom=0.05e5 "nominal pressure drop";
    parameter SI.VolumeFlowRate Vdot_nom=0.4 "nominal volume flow rate";
    final parameter Real K=(deltaP_nom/Vdot_nom)/n;

  equation
    // die max-Abfrage braucht man, für den Fall, daß man die Kolonne leer initialisiert und das Inertgas nicht mit abbildet: in dem Fall würde sonst ein Volumenstrom in die falsche Richtung entstehen

    for j in 1:n - 1 loop
        Vdot[j] = if startUp[j] then 0 else  max(0,sign(p[j] - p[j + 1])*max(0, (abs(p[j] - p[j + 1])/K)));
    end for;

    /*** outlet: half of the n-th discreet element ***/
    Vdot[n] = if startUp[n] then 0 else max(0,sign(p[n] - p[n + 1])*max(0, (abs(p[n] - p[n + 1])/(K/2))));

    /*** entry: half of the first discreet element ***/
    p_v_in = Vdot_in*K/2 + p[1];

  assert(not homotopyMethod.bool_dp,"this pressure loss model does not support homotopy - please deactivate");

    annotation (Documentation(info="<html>
<p>Linear pressure drop, using a nominal value for delta_P and Vdot </p>
</html>"));
  end NominalLinear2;
end TrayColumn;
