within ThermalSeparation.FilmModel.StructuredPackedColumn;
partial model BaseFilmPacked "base film model for packed column"
  replaceable record Geometry =
      ThermalSeparation.Geometry.StructuredPackedColumn.Mellapak250Y                constrainedby
    ThermalSeparation.Geometry.StructuredPackedColumn.Geometry
      annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

end BaseFilmPacked;
