within ThermalSeparation.HeatAndMassTransfer.TrayColumn.Vapour;
partial model BaseVapMT "base model for vapour mass transfer coefficient"
  parameter Integer n(min=1) annotation(Dialog(enable=false));
  parameter Integer n_k
    "number of calculated mass transfer coefficients - nSL or k";
        replaceable record Geometry =
      ThermalSeparation.Geometry.PlateColumn.Geometry constrainedby
    ThermalSeparation.Geometry.PlateColumn.Geometry;
    Geometry geometry(final n=n);
    input SI.Height h[n] "height of the two-phase regime on the tray";
input Real eps_liq_2ph[n]
    "fraction of liquid in the two-phase regime on the tray";
    input SI.DiffusionCoefficient D[ n,n_k];
  input SI.Velocity w_sup[n] "superficial vapour velocity";
  output ThermalSeparation.Units.CoefficentOfMassTransfer k[n,n_k];
equation

end BaseVapMT;
