within ThermalSeparation.Holdup;
package TrayColumn
  partial model BaseHoldup
    parameter Integer n(min=1) annotation(Dialog(enable = false));
    parameter Boolean considerStartUp annotation(Dialog(enable = false));
    parameter Boolean constProp = false
      "constant values for density and surface tension";
    parameter SI.Density rho_l_const = 1000 "constant value for liquid density" annotation(Dialog(enable = constProp));
    parameter SI.Density rho_v_const = 2 "constant value for vapour density" annotation(Dialog(enable = constProp));
    parameter SI.SurfaceTension sigma_const = 0.075
      "constant value for surface tension"                                               annotation(Dialog(enable = constProp));
    input SI.Density rho_l[n];
    input SI.Density rho_v[n];
    input SI.SurfaceTension sigma[n];
    input Real omega[n];
    input SI.Length h[n] "height of the 2ph region on the tray";
    input SI.VolumeFlowRate Vdot_v[n] "vapour volume flow rate";
    output SI.VolumeFlowRate Vdot[n];
    output Real eps_liq_2ph[n](start=0.5*ones(n))
      "liquid fraction in the two-phase area on the plate";
    output ThermalSeparation.Units.F_Factor F[n](
                           start=ones(n)*200)
      "F-factor based on the active area";
    output ThermalSeparation.Units.F_Factor F_max[n] "maximum vapour load";

    replaceable record Geometry =
        ThermalSeparation.Geometry.PlateColumn.Geometry;
    Geometry geometry(n=n);

    SI.Density rho_liq[n] = if constProp then fill(rho_l_const,n) else rho_l;
    SI.Density rho_vap[n] = if constProp then fill(rho_v_const,n) else rho_v;
    SI.SurfaceTension sigma_liq[n] = if constProp then fill(sigma_const,n) else sigma;
  equation

  end BaseHoldup;

  model Stichlmair "Stichlmair"
    extends BaseHoldup;

      //minimum and maximum vapour load
    parameter ThermalSeparation.Units.F_Factor F_max_const = 1
      "guess initial max. vapour load";
     parameter ThermalSeparation.Units.F_Factor F_const = 0.2
      "guess initial vapour load";
     parameter Real time_F = 0
      "time, when the F-factor shall be no longer constant, but shall be calculated";

  protected
     Real aux_F = min(1,max(0,0.1*time - time_F))
      "auxilary variable to change the F-factor from a constant value to a value calculated by a correlation";
    ThermalSeparation.Units.F_Factor F_corr[n]
      "correlation for F-factor based on the active area";
    ThermalSeparation.Units.F_Factor F_max_corr[n]
      "correlation for F_max,  eq. from Stichlmair, Chem.Ing.Tech. 50 (1978) Nr. 4";
    ThermalSeparation.Units.F_Factor F_min[n];
    ThermalSeparation.Units.F_Factor F_min_1[n]
      "minimum vapour load, eq. from Mersmann, Chem.Ing.Tech. 35 (1963) Nr. 2";
    ThermalSeparation.Units.F_Factor F_min_2[n]
      "minimum vapour load, eq. from Ruff, Chem.Ing.Techn. 48 (1976) Nr.9";

   SI.Length delta_h_corr[n]
      "correction term for the calculation of the height of the two-phase regime, Stichlmair";

  equation
  for j in 1:n loop
  /*** during startUp there is no, or nearly no vapour in the liquid on the stage, since all vapour is condensing directly ***/
      eps_liq_2ph[j] = 0.99*(1-omega[j])+  max(1e-5,1 - (F[j]/F_max[j])^0.28)*omega[j];
      F_max[j] =  F_max_corr[j]*aux_F + F_max_const*(1-aux_F);
      F_max_corr[j] =  2.5 * ((geometry.phi)^2 * sigma_liq[j]*(rho_liq[j] - rho_vap[j])*Modelica.Constants.g_n)^(1/4);
      F[j] = F_corr[j]*aux_F + F_const *(1-aux_F);
      F_corr[j] =  max(1e-5,Vdot_v[j]/geometry.A_plate * sqrt(rho_vap[j]));
      F_min_1[j] = sqrt(2*sigma_liq[j]/geometry.d_h);
      F_min_2[j] = sqrt(0.37 * geometry.d_h * Modelica.Constants.g_n * (rho_liq[j] - rho_vap[j])^(5/4)/rho_vap[j]^0.25);
      F_min[j] = min(F_min_1[j],F_min_2[j]);

        /*** Frickel-Lsg, da sonst die StartUp-Rektifikationskolonne nicht läuft ***/
      delta_h_corr[j] =if considerStartUp then 0 else omega[j]* max(0,125/(Modelica.Constants.g_n * (rho_liq[j]-rho_vap[j])) * (F[j] - 0.2 * sqrt(rho_vap[j]) / (1-eps_liq_2ph[j]))^2);

  /*** Stichlmair: Dimensionierung des Gas/Flüssigkeits-Kontaktapparates Bodenkolonne, Teil II, Chem. Ing. Tech. 50, Nr. 5, Gl. (3.2) ***/
     Vdot[j] = if h[j]/eps_liq_2ph[j] > geometry.h_w then sign(h[j]/eps_liq_2ph[j]-geometry.h_w)* geometry.l_w
   *eps_liq_2ph[j]*(abs(h[j]/eps_liq_2ph[j]-geometry.h_w - delta_h_corr[j])*Modelica.Constants.g_n^(1/3)
   /1.45)^(3/2) else 0;
  end for;

    annotation (Documentation(info="<html>
<p>Correlation for liquid volume flow rate: equation (3.2) from publication [1].</p>
<pre> 
<font style=\"color: #006400; \">[1] Stichlmair:&nbsp;Dimensionierung&nbsp;des&nbsp;Gas/Fl&uuml;ssigkeits-Kontaktapparates&nbsp;Bodenkolonne,&nbsp;Teil&nbsp;II,&nbsp;Chem.&nbsp;Ing.&nbsp;Tech.&nbsp;50,&nbsp;Nr.&nbsp;5</font></pre>
</html>"));
  end Stichlmair;

  model Stichlmair_simplified
    "Stichlmair simplified - without vapour in liquid regime"
    extends BaseHoldup;

      //minimum and maximum vapour load
    parameter ThermalSeparation.Units.F_Factor F_max_const = 1
      "guess initial max. vapour load";
     parameter ThermalSeparation.Units.F_Factor F_const = 0.2
      "guess initial vapour load";
     parameter Real time_F = 0
      "time, when the F-factor shall be no longer constant, but shall be calculated";

  protected
     Real aux_F = min(1,max(0,0.1*time - time_F))
      "auxilary variable to change the F-factor from a constant value to a value calculated by a correlation";
    ThermalSeparation.Units.F_Factor F_corr[n]
      "correlation for F-factor based on the active area";
    ThermalSeparation.Units.F_Factor F_max_corr[n]
      "correlation for F_max,  eq. from Stichlmair, Chem.Ing.Tech. 50 (1978) Nr. 4";

   SI.Length delta_h_corr[n]
      "correction term for the calculation of the height of the two-phase regime, Stichlmair";

  equation
  for j in 1:n loop
  /*** during startUp there is no, or nearly no vapour in the liquid on the stage, since all vapour is condensing directly ***/
      eps_liq_2ph[j] =1;// 0.99*(1-omega[j])+  max(1e-5,1 - (F[j]/F_max[j])^0.28)*omega[j];
      F_max[j] =  F_max_corr[j]*aux_F + F_max_const*(1-aux_F);
      F_max_corr[j] =  2.5 * ((geometry.phi)^2 * sigma_liq[j]*(rho_liq[j] - rho_vap[j])*Modelica.Constants.g_n)^(1/4);
      F[j] = F_corr[j]*aux_F + F_const *(1-aux_F);
      F_corr[j] =  max(1e-5,Vdot_v[j]/geometry.A_plate * sqrt(rho_vap[j]));

        /*** Frickel-Lsg, da sonst die StartUp-Rektifikationskolonne nicht läuft ***/
      delta_h_corr[j] =if considerStartUp then 0 else omega[j]* max(0,125/(Modelica.Constants.g_n * (rho_liq[j]-rho_vap[j])) * (F[j] - 0.2 * sqrt(rho_vap[j]) / (1-eps_liq_2ph[j]))^2);

  /*** Stichlmair: Dimensionierung des Gas/Flüssigkeits-Kontaktapparates Bodenkolonne, Teil II, Chem. Ing. Tech. 50, Nr. 5, Gl. (3.2) ***/
     Vdot[j] = if h[j]/eps_liq_2ph[j] > geometry.h_w then sign(h[j]/eps_liq_2ph[j]-geometry.h_w)* geometry.l_w
   *eps_liq_2ph[j]*(abs(h[j]/eps_liq_2ph[j]-geometry.h_w - delta_h_corr[j])*Modelica.Constants.g_n^(1/3)
   /1.45)^(3/2) else 0;
  end for;

    annotation (Documentation(info="<html>
<p>Correlation for liquid volume flow rate: equation (3.2) from publication [1].</p>
<pre> 
<font style=\"color: #006400; \">[1] Stichlmair:&nbsp;Dimensionierung&nbsp;des&nbsp;Gas/Fl&uuml;ssigkeits-Kontaktapparates&nbsp;Bodenkolonne,&nbsp;Teil&nbsp;II,&nbsp;Chem.&nbsp;Ing.&nbsp;Tech.&nbsp;50,&nbsp;Nr.&nbsp;5</font></pre>
</html>"));
  end Stichlmair_simplified;
end TrayColumn;
