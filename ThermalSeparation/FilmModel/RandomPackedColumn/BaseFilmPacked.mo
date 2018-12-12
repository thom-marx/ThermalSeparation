within ThermalSeparation.FilmModel.RandomPackedColumn;
partial model BaseFilmPacked "base film model for packed column"
  replaceable record Geometry =
      ThermalSeparation.Geometry.RandomPackedColumn.Geometry                constrainedby
    ThermalSeparation.Geometry.RandomPackedColumn.Geometry
      annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

end BaseFilmPacked;
