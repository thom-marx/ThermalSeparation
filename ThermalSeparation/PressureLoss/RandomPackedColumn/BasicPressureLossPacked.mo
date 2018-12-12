within ThermalSeparation.PressureLoss.RandomPackedColumn;
partial model BasicPressureLossPacked
  "basic pressure loss model for a plate column"
  extends BasicPressureLoss;
    replaceable record Geometry =
      ThermalSeparation.Geometry.RandomPackedColumn.Geometry                            constrainedby
    ThermalSeparation.Geometry.RandomPackedColumn.Geometry annotation(Dialog(enable=false));
    Geometry geometry(n=n);
  input Real hu_dyn[n] "dynamic holdup";

end BasicPressureLossPacked;
