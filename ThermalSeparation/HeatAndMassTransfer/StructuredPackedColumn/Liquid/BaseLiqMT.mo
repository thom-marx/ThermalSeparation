within ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Liquid;
partial model BaseLiqMT "base model for liquid mass transfer coefficient"
  parameter Integer n(min=1) annotation(Dialog(enable=false));
  parameter Integer n_k(min=1)
    "number of calculated mass transfer coefficients - nSL or k";

          replaceable record Geometry = ThermalSeparation.Geometry.StructuredPackedColumn.Geometry
                                                                 constrainedby ThermalSeparation.Geometry.StructuredPackedColumn.Geometry;
    Geometry geometry(final n=n);
      input SI.SurfaceTension sigma[n];
      input SI.Velocity w_sup[n] "superficial vapour velocity";
     input SI.DynamicViscosity eta[n];
  input SI.Density rho[n];
  input SI.DiffusionCoefficient D[n,n_k];
  input Real eps_liq[n];
  output ThermalSeparation.Units.CoefficentOfMassTransfer k[n,n_k];
equation

end BaseLiqMT;
