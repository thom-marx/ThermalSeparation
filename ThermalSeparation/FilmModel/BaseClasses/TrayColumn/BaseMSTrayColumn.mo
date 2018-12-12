within ThermalSeparation.FilmModel.BaseClasses.TrayColumn;
partial model BaseMSTrayColumn
  "base model for Maxwell-Stefan R-matrix in a tray column"
  extends ThermalSeparation.FilmModel.BaseClasses.ComputeR;

replaceable record Geometry = 
      ThermalSeparation.Geometry.PlateColumn.Geometry                            constrainedby
    ThermalSeparation.Geometry.PlateColumn.Geometry;
    Geometry geometry(final n=n);

end BaseMSTrayColumn;
