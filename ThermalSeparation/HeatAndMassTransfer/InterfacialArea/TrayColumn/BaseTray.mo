within ThermalSeparation.HeatAndMassTransfer.InterfacialArea.TrayColumn;
model BaseTray
  extends Base;

  replaceable record Geometry =
      ThermalSeparation.Geometry.PlateColumn.Geometry                            constrainedby ThermalSeparation.Geometry.PlateColumn.Geometry;
    Geometry geometry(final n=n);
  input Real F[n];
  input Real F_max[n];
  input SI.Density rho_v[n] "vapour density";

end BaseTray;
