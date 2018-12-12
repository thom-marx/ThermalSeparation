within ThermalSeparation.FilmModel.BaseClasses.RandomPackedColumn;
partial model BaseMSPackedColumn
  "base model for Mawell-Stefan R-matrix in a packed column"
  extends ThermalSeparation.FilmModel.BaseClasses.ComputeR;

input Real eps_liq[n];

      replaceable record Geometry = 
      ThermalSeparation.Geometry.RandomPackedColumn.Geometry                  constrainedby
    ThermalSeparation.Geometry.RandomPackedColumn.Geometry;
    Geometry geometry(final n=n);
equation

end BaseMSPackedColumn;
