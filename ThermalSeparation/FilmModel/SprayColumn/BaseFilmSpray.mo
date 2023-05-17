within ThermalSeparation.FilmModel.SprayColumn;
partial model BaseFilmSpray "base film model for spray column"
  parameter Integer nn(min=1)   annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

   input SI.Diameter d_drop[nn](start=1e-3*ones(nn));
   input Real n_drop[nn](start=1*ones(nn));
       replaceable record Geometry = ThermalSeparation.Geometry.SprayColumn.GeometrySpray
                                                           constrainedby ThermalSeparation.Geometry.SprayColumn.Geometry
            annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

  Geometry geometry(n=nn);
equation

end BaseFilmSpray;
