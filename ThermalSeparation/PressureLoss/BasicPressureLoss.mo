within ThermalSeparation.PressureLoss;
partial model BasicPressureLoss

  replaceable package MediumLiquid = Media.BaseMediumLiquid;
  replaceable package MediumVapour = Media.BaseMediumVapour;
  input MediumLiquid.ThermodynamicProperties propsLiq[n];
  input MediumVapour.ThermodynamicProperties propsVap[n];
  input SI.VolumeFlowRate Vdot_l[n];
  parameter Integer n(min=1)=8;
  input SI.Pressure p[n+1];
  input SI.VolumeFlowRate Vdot_in;
  input Real eps_liq[n];

  parameter Boolean rho_const = false
    "rho is set constant (rho = rho_nom) for calculation of k" annotation(Dialog(group="Use of nominal values"));
  parameter SI.Density rho_l_nom = 1000
                                 annotation(Dialog(enable=rho_const,group="Use of nominal values"));
  parameter SI.Density rho_v_nom = 1
                                 annotation(Dialog(enable=rho_const,group="Use of nominal values"));

  output SI.VolumeFlowRate Vdot[n];
  output SI.Pressure p_v_in;
//  Boolean useHomotopy;
  replaceable model HomotopyMethod =
      ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.NoHomotopy
                       constrainedby ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.BaseHomotopy
                                                                                      annotation(Dialog(enable=false,choicesAllMatching=true));

HomotopyMethod homotopyMethod(n=n);

protected
  SI.Density rho_l[n] = propsLiq.rho;
  SI.SurfaceTension sigma_l[n] = propsLiq.sigma;
  SI.Density rho_v[n] = propsVap.rho;
  SI.Density rho_l_calc[n] = if rho_const then fill(rho_l_nom,n) else rho_l;
  SI.Density rho_v_calc[n] = if rho_const then fill(rho_v_nom,n) else rho_v;

end BasicPressureLoss;
