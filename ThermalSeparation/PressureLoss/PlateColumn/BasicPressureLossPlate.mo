within ThermalSeparation.PressureLoss.PlateColumn;
partial model BasicPressureLossPlate
  "basic pressure loss model for a plate column"
  extends BasicPressureLoss;
    replaceable record Geometry = 
      ThermalSeparation.Geometry.PlateColumn.Geometry                            constrainedby
    ThermalSeparation.Geometry.PlateColumn.Geometry annotation(Dialog(enable=false));
    Geometry geometry(n=n);

  input SI.Height h[n](start=0.9*ones(n))
    "height of the 2ph region on the tray";
  input Real eps_liq_2ph[n] "liquid fraction in the two-phase area on the tray";
equation

end BasicPressureLossPlate;
