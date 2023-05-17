within ThermalSeparation.FilmModel.SprayColumn;
model EquilibriumEffectiveDiff "equilibrium on each stage 3"
  import ThermalSeparation;
  extends BaseNonEqType(nn=n); //for replaceability
  extends ThermalSeparation.FilmModel.BaseClasses.EquilibriumEffectiveDiff(
                                                               redeclare record BaseGeometry =
        Geometry);
        Geometry geometry;
equation
    /***for equilibrium model of spray column: deviation from equilibrium calculated using HETP values***/
A_I=fill(400*geometry.H*geometry.A/n,n);
  for j in 1:n loop
    for i in 1:nSV loop
       x_v_star[j,i] =  x_v_star_eq[j,i];
    end for;
    end for;
end EquilibriumEffectiveDiff;
