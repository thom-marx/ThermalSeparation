within ThermalSeparation.HeatAndMassTransfer.PlateColumn;
partial model Base
  extends ThermalSeparation.HeatAndMassTransfer.ComputeR;

replaceable record Geometry = 
      ThermalSeparation.Geometry.PlateColumn.Geometry                            constrainedby
    ThermalSeparation.Geometry.PlateColumn.Geometry;
    Geometry geometry(final n=n);

end Base;
