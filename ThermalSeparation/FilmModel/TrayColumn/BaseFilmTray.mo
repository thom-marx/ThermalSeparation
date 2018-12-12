within ThermalSeparation.FilmModel.TrayColumn;
partial model BaseFilmTray "base film model for tray column"
    parameter Integer nn(min=1)
     annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

  replaceable record Geometry =
      ThermalSeparation.Geometry.PlateColumn.Geometry                           constrainedby
    ThermalSeparation.Geometry.PlateColumn.Geometry
  annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

  input Real eps_liq_2ph[nn]
    "liquid fraction in the two-phase area on the plate";
  input ThermalSeparation.Units.F_Factor F[nn]
    "F-factor based on the active area";
  input ThermalSeparation.Units.F_Factor F_max[nn] "maximum vapour load";
  input SI.Height h[nn];
end BaseFilmTray;
