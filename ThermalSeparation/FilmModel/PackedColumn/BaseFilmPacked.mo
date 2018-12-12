within ThermalSeparation.FilmModel.PackedColumn;
partial model BaseFilmPacked "base film model for packed column"
  replaceable record Geometry =
      ThermalSeparation.Geometry.PackedColumn.Geometry                          constrainedby
    ThermalSeparation.Geometry.PackedColumn.Geometry
      annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

end BaseFilmPacked;
