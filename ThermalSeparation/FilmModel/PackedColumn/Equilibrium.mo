within ThermalSeparation.FilmModel.PackedColumn;
model Equilibrium "equilibrium on each stage"

  extends BaseFilmPacked; //for replaceability
  extends ThermalSeparation.FilmModel.BaseClasses.Equilibrium;
equation
    /***for equilibrium model of packed column: deviation from equilibrium calculated using HETP values***/

  for j in 1:n loop
    for i in 1:nSV loop
       x_v_star[j,i] =  x_v_star_eq[j,i];
    end for;
    end for;
end Equilibrium;
