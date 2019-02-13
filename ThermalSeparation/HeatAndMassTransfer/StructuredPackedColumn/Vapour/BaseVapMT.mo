within ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Vapour;
partial model BaseVapMT "base model for vapour mass transfer coefficient"
  parameter Integer n(min=1) annotation(Dialog(enable=false));
  parameter Integer n_k
    "number of calculated mass transfer coefficients - nSL or k";
              replaceable record Geometry =
      ThermalSeparation.Geometry.StructuredPackedColumn.Geometry constrainedby
    ThermalSeparation.Geometry.StructuredPackedColumn.Geometry;
    Geometry geometry(final n=n);
      input SI.SurfaceTension sigma[n];
         input SI.DynamicViscosity eta[n];
  input SI.Density rho[n];
  input SI.DiffusionCoefficient D[n,n_k];
    input Real eps_liq[n];
     input SI.Velocity w_sup_l[   n] "superficial liquid velocity";
    input SI.Velocity w_sup_v[n] "superficial vapour velocity";
  output ThermalSeparation.Units.CoefficentOfMassTransfer k[n,n_k];
equation

end BaseVapMT;
