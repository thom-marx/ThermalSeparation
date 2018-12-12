within ThermalSeparation;
package FilmModel 
annotation (Documentation(info="<html>
<p>The UML-diagram shows the structure of the film models. Different film models are located in the package <a href=\"modelica://ThermalSeparation.FilmModel.BaseClasses\">BaseClasses</a>, namely film models for <a href=\"modelica://ThermalSeparation.FilmModel.BaseClasses.MaxwellStefan\">Maxwell-Stefan Mass Transfer</a>, <a href=\"modelica://ThermalSeparation.FilmModel.BaseClasses.TrueEquilibrium\">True equilibrium</a>, <a href=\"modelica://ThermalSeparation.FilmModel.BaseClasses.Enhancement\">Enhancement factor</a>, ... </p>
<p>Most film models require values for variables such as interfacial area, mass or heat transfer coefficient. Correlations for such variables are often column type, depended, therefore there exist a package for each column type where the base film models are extended (boxes marked in grey). These models supply the equations for column specific variables. A common base class such as <a href=\"modelica://ThermalSeparation.FilmModel.StructuredPackedColumn.BaseFilmPacked\">BaseFilmSP</a> (marked with grey hatches) groups all film models referring to one column type.</p>
<p><br/><img src=\"modelica://ThermalSeparation/Pictures/UML_BaseFilm2.png\"/></p>
</html>"));
end FilmModel;
