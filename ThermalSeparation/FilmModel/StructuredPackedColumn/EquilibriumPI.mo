within ThermalSeparation.FilmModel.StructuredPackedColumn;
model EquilibriumPI "equilibrium on each stage using a PI controller"
  import ThermalSeparation;

  extends BaseNonEqType; //for replaceability
  extends ThermalSeparation.FilmModel.BaseClasses.EquilibriumPI(
                                                              redeclare record BaseGeometry =
        ThermalSeparation.Geometry);
        Geometry geometry;
equation
   A_I = fill(geometry.a * geometry.H*geometry.A /n,n);

    /***for equilibrium model of packed column: deviation from equilibrium calculated using HETP values***/

  for j in 1:n loop
    for i in 1:nSV loop
       x_v_star[j,i] =  x_v_star_eq[j,i];
    end for;
    end for;
end EquilibriumPI;
