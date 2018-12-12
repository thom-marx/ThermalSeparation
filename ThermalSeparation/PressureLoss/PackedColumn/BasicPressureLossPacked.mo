within ThermalSeparation.PressureLoss.PackedColumn;
partial model BasicPressureLossPacked
  "basic pressure loss model for a plate column"
  extends BasicPressureLoss;
    replaceable record Geometry = 
      ThermalSeparation.Geometry.PackedColumn.Geometry                            constrainedby
    ThermalSeparation.Geometry.PackedColumn.Geometry annotation(Dialog(enable=false));
    Geometry geometry(n=n);
  input Real hu_dyn[n] "dynamic holdup";

end BasicPressureLossPacked;
