within ThermalSeparation.FilmModel.SprayColumn;
model TrueEquilibrium "equilibrium on each stage 2"
  import ThermalSeparation;
  extends BaseEquilibriumType(nn=n); //for replaceability
  extends ThermalSeparation.FilmModel.BaseClasses.TrueEquilibrium(
                                                               redeclare record
      BaseGeometry =
        ThermalSeparation.Geometry);
equation
    /***for equilibrium model of spray column: deviation from equilibrium calculated using HETP values***/

  for j in 1:n loop
    for i in 1:nSV loop
       x_v_star[j,i] =  x_v_star_eq[j,i];
    end for;
    end for;
end TrueEquilibrium;
