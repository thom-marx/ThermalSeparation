within ThermalSeparation.FilmModel.BaseClasses.PackedColumn;
partial model BaseMSPackedColumn
  "base model for Mawell-Stefan R-matrix in a packed column"
  extends ThermalSeparation.FilmModel.BaseClasses.ComputeR;

input Real eps_liq[n];

      replaceable record Geometry = 
      ThermalSeparation.Geometry.StructuredPackedColumn.Geometry                  constrainedby
    ThermalSeparation.Geometry.StructuredPackedColumn.Geometry;
    Geometry geometry(final n=n);
equation

end BaseMSPackedColumn;
