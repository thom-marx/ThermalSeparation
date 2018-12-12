within ThermalSeparation.FilmModel.StructuredPackedColumn;
model EquilibriumBalancing "EquilibriumBalancing"
  import ThermalSeparation;

  extends BaseEquilibriumType; //for replaceability
  extends ThermalSeparation.FilmModel.BaseClasses.EquilibriumBalancing(
       redeclare record BaseGeometry =           ThermalSeparation.Geometry);
equation
    /***for equilibrium model of packed column: deviation from equilibrium calculated using HETP values***/

  for j in 1:n loop
    for i in 1:nSV loop
       x_v_star[j,i] =  x_v_star_eq[j,i];
    end for;
    end for;
end EquilibriumBalancing;
