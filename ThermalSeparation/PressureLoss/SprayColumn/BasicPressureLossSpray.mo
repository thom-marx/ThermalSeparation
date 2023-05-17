within ThermalSeparation.PressureLoss.SprayColumn;
partial model BasicPressureLossSpray
  "basic pressure loss model for a spray column"
  extends BasicPressureLoss;
    replaceable record Geometry =
      ThermalSeparation.Geometry.SprayColumn.Geometry                            constrainedby ThermalSeparation.Geometry.SprayColumn.Geometry
                                                     annotation(Dialog(enable=false));
    Geometry geometry(n=n);

end BasicPressureLossSpray;
