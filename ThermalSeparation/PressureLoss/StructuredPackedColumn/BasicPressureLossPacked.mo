within ThermalSeparation.PressureLoss.StructuredPackedColumn;
partial model BasicPressureLossPacked
  "basic pressure loss model for a plate column"
  extends BasicPressureLoss;
    replaceable record Geometry =
      ThermalSeparation.Geometry.StructuredPackedColumn.Geometry                  constrainedby
    ThermalSeparation.Geometry.StructuredPackedColumn.Geometry
                                                     annotation(Dialog(enable=false));
    Geometry geometry(n=n);
  input Real hu_dyn[n] "dynamic holdup";

end BasicPressureLossPacked;
