within ThermalSeparation;
package UsersGuide "Users Guide"


annotation(DocumentationClass=true, Documentation(info="<html>
<p><u><b><span style=\"color: #008000;\">Solver Type</span></b></u></p>
<p>In general Dassl works best.</p>
<p><u><b><span style=\"color: #008000;\">Test Models</span></b></u></p>
<p>The functionality of this library was tested with models which are located in the package <a href=\"Modelica://ThermalSeparation.Examples\">Examples</a>. It&apos;s a good start to study this models to gain experience in using this library.</p>
<p><u><b><span style=\"color: #008000;\">Column Types</span></b></u></p>
<p>Four different column types can be found in the package <a href=\"Modelica://ThermalSeparation.Components.Columns\">Components.Columns</a>, namely tray column, spray column, random packed column and structured packed column.</p>
<p><u><b><span style=\"color: #008000;\">Feed</span></b></u></p>
<p>There are no special feed segments, each discrete element could be in theory a feed tray. If a column has a feed, the boolean parameter <i><b>hasLiquidFeed</b></i> and/or <i><b>hasVapourFeed</b></i> in the tab <i>Feed</i> has to be set to true (default value: false). The user has then to supply the number of feeds the column has (<i><b>numberLiquidFeeds</b></i> / <i><b>numberVapourFeeds</b></i>). And the number of the stage, where the feed enters the column (<i><b>stageLiquidFeed</b></i> / <i><b>stageVapourFeed</b></i>). The numbering of stages where the feeds enter the column are provided as a vector (like <i>stageLiquidFeed</i> = {3,5}), this makes it possible to have more than one feed. The stage at the bottom of the column is the first stage.</p>
<p><u><b><span style=\"color: #008000;\">Medium Models </span></b></u></p>
<p>Vapour and liquid medium packages have to be provided separately.</p>
<p>The liquid medium models and the vapour medium models can differ both in the number of mediums they contain as well as in the substance types. The parameter <i><b>nS </b></i>is the number of substances which are in the liquid as well as in the vapour phase. This value has to be supplied on the top-level in the GUI of the columns. The arrangement of the different substances in the medium models is arbitrary. The parameter <i><b>mapping </b></i>has to be used to map the different vectors one to another. </p>
<p><u><b>Example:</b></u> Vapour = {N2, H2O, CO2}, Liquid = {N2, H+, HCO3-, H2O, CO2} , mapping = {{1,1},{2,4},{3,5}}.</p>
<p>The value for the parameter<i><b> mapping</b></i> is also supplied on the top-level of the column GUI.</p>
<p>In the vector <i><b>inertVapour </b></i>and<i><b> inertLiquid </b></i>it has to be stated whether the substance is to be considered as inert substance or not. For the example above:</p>
<p>inertVapour = {false, false, false}, inertLiquid = {false, true, true, false, false}.</p>
<p>If new medium models are added to the library they have to extend from <a href=\"Modelica://ThermalSeparation.Media.BaseMediumVapour\">BaseMediumVapour</a> and <a href=\"Modelica://ThermalSeparation.Media.BaseMediumLiquid\">BaseMediumLiquid</a> respectively (or <a href=\"Modelica://ThermalSeparation.Media.BaseMediumLiquidReaction\">BaseMediumLiquidReactions</a>, if reaction will occur in the liquid phase).</p>
<p><u><b><span style=\"color: #008000;\">Film Models</span></b></u></p>
<p>The film model provides the relation between the bulk mole fractions and the molar flow rates between the bulk phases and the film phases.</p>
<p>For a detailed description of the film model classes, see package <a href=\"modelica://ThermalSeparation.FilmModel\">FilmModel.</a></p>
<p>The film model is chosen on top-level of the column GUI. Most film models require values for variables such as interfacial area, mass or heat transfer coefficient. Correlations for such variables are not determined on top level in the GUI, but when clicking on the modifier of the film model. In this case the redeclaration window will open and correlations for heat and mass transfer and so on can be chosen (if applicable for the chosen film model). This can be seen in the following screenshot:</p>
<p><br><img src=\"modelica://ThermalSeparation/Pictures/DymolaScreenshot5.jpg\"/></p>
<p>In the film model also the numeric states can be chosen (tab &QUOT;General&QUOT; of the redeclaration window). Different suggestions of state variables are selected in the package <a href=\"modelica://ThermalSeparation.FilmModel.StateSelection\">StateSelection</a> where also own state options can be added.</p>
<p><u><b><span style=\"color: #008000;\">Initialization</span></b></u></p>
<p>The Initialization in steady-state is not (yet) possible, and different sets of initial equations exist, hoping that one of them will work. The set of initial equations used can be selected in the tab &QUOT;Initializiation&QUOT; of the GUI. Several variables are listed in this tab, for which start values are required. This does not necessarily mean that these values will be physical start values of t=0s. They may also be only iteration start values. Which variables are physical start values and which variables are separation start values depend on the initial equations chosen. Own initial equations can be added and have to extend from the model <a href=\"modelica://ThermalSeparation.Components.Columns.BaseClasses.Initialization.BaseInit\">BaseInit</a>.</p>
</html>"));
end UsersGuide;
