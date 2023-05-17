within ThermalSeparation.FilmModel;
package StructuredPackedColumn
  model TrueEquilibriumStartUpCCSAbsorption "true equilibrium start up for absorption process"

    extends BaseFilmPacked; //for replaceability
    extends ThermalSeparation.FilmModel.BaseClasses.TrueEquilibriumStartUpCCSAbsorption(
         redeclare record BaseGeometry = Geometry);
  equation
      /***for equilibrium model of packed column: deviation from equilibrium calculated using HETP values***/

    for j in 1:n loop
      for i in 1:nSV loop
         x_v_star[j,i] =  x_v_star_eq[j,i];
      end for;
      end for;
  end TrueEquilibriumStartUpCCSAbsorption;

  model TrueEquilibriumStartUpCCSDesorption "true equilibrium start up for desorption process"

    extends BaseFilmPacked; //for replaceability
    extends ThermalSeparation.FilmModel.BaseClasses.TrueEquilibriumStartUpCCSDesorption(
         redeclare record BaseGeometry = Geometry);
  equation
      /***for equilibrium model of packed column: deviation from equilibrium calculated using HETP values***/

    for j in 1:n loop
      for i in 1:nSV loop
         x_v_star[j,i] =  x_v_star_eq[j,i];
      end for;
      end for;
  end TrueEquilibriumStartUpCCSDesorption;
end StructuredPackedColumn;
