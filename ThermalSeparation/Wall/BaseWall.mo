within ThermalSeparation.Wall;
partial model BaseWall "Model of the heat loss at the column outer wall"
  parameter Integer n(min=1) annotation(Dialog(enable=false));
  parameter Boolean stat=false "true if wall is calculated stationary";
  input SI.Temperature T_amb;

  replaceable record Geometry =
      ThermalSeparation.Geometry.BasicGeometry                          constrainedby ThermalSeparation.Geometry.BasicGeometry  annotation(choicesAllMatching);
  Geometry geometry(n=n) "Geometry Record of the corresponding section";
  SI.Temperature T_wall[n] "wall temperature";
  parameter SI.Temperature T_start "initial temperature of the wall" annotation(Dialog(enable=false));
  input SI.Temperature T[n];
  output SI.HeatFlowRate Qdot[n]
    "heat flow from inner side of the wall to the wall";

  /*** to be supplied in extending class: ***/
  SI.HeatFlowRate Qdot_out[n]
    "heat flow from the wall to the outer side of the wall";
  parameter SI.SurfaceCoefficientOfHeatTransfer alpha_inner= 200
    "heat transfer coefficient at inner side of the wall";
//Dummy temperature which is needed in order to use the temperature as state if the wall is not calculated steady-state
protected
  SI.Temperature T_dummy[n](each stateSelect=StateSelect.prefer, start=fill(T_start,n)) = T_wall if not stat;

equation
  for j in 1:n loop
    Qdot[j]=alpha_inner * 2 * Modelica.Constants.pi * geometry.H /n * geometry.r1 * (T[j]-T_wall[j]);
  if stat then
    0 = Qdot[j] - Qdot_out[j];
  else
    geometry.rho_wall*geometry.H/n*Modelica.Constants.pi*(geometry.r2^2 - geometry.r1^2)*geometry.c_wall*der(T_wall[j]) = Qdot[j] - Qdot_out[j];
  end if;
  end for;

initial equation
    if not stat then
    T_wall = fill(T_start,n);
    end if;

  annotation (Diagram(graphics),
                       Icon(graphics={
        Rectangle(
          extent={{-40,80},{40,-80}},
          lineColor={0,0,0},
          lineThickness=1,
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Text(
          extent={{42,50},{98,-48}},
          lineColor={0,0,0},
          lineThickness=1,
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward,
          textString=
               "Tamb"),
        Line(
          points={{-40,52},{-32,14},{-26,2},{-16,-10},{-2,-22},{10,-28},{26,-34},
              {34,-36},{40,-36}},
          color={255,0,0},
          thickness=1)}),
    Documentation(revisions="<html>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"1\"><tr>
<td><p align=\"center\"><h4>created by </h4></p></td>
<td><p><a href=\"mailto:karin.dietl@tu-harburg.de\">Karin Dietl</a> &AMP; <a href=\"mailto:andreas.joos@tu-harburg.de\">Andreas Joos</a></p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>creation date </h4></p></td>
<td><p>01.01.2009 </p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>revised by </h4></p></td>
<td><p>nobody so far</p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>last revision </h4></p></td>
<td><p>this is an alpha version... </p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>based on </h4></p></td>
<td></td>
</tr>
</table>
<p><br/><br/><br/>Documentation last revised: 20.7.2011</p>
</html>", info="<html>
<p>This is the partial model of the column wall. The heat resistance in the column wall is neglected. The heat resistance between column wall and liquid phase and between column wall and ambiance is taken into account. If there exist an insulation, the heat resistance through the insulation is also taken into account. There exists the possiblity to neglect the heat capacitance of the wall, in this case the boolean parameter <code>stat</code> must be set to true. The initial temperature of the wall is supposed to be equal to the initial temperature of the liquid in the column. The insulation is always calculated steady-state, since it is supposed that the heat capacitance of the insulation is always much smaller than the heat capacitance of the column wall.</p>
<p>The heat resistance at the column inner wall can only be taken into account using a constant value for the heat transfer coefficient alpha. </p>
<p><br/>The extending class has to supply an equation for the heat flow rate between insulation (or column wall, if the column has no insulation at all) and ambiance. In the class <a href=\"Modelica://ThermalSeparation.Wall.Adiabatic\">Adiabatic</a> it is assumed that this heat flow rate is zero. <a href=\"Modelica://ThermalSeparation.Wall.ConstAlpha\">ConstAlpha</a> assumes a constant heat transfer coeffient to calculate the heat flow rate and <a href=\"Modelica://ThermalSeparation.Wall.NaturalConvection\">NaturalConvection</a> calculates the heat transfer coefficent using correlations for natural convection. <a href=\"Modelica://ThermalSeparation.Wall.ConstQdot\">Qdot</a> assumes a constant flow rate to ambiance.</p>
<p>Technische Universit&auml;t Hamburg-Harburg </p>
<p>Institut f&uuml;r Thermofluiddynamik, Technische Thermodynamik </p>
<p>Denickestra&szlig;e 17 </p>
<p>D-21073 Hamburg </p>
<p>Germany </p>
<p><b><a href=\"http://www.tu-harburg.de/tt\">www.tu-harburg.de/tt</a></b></p>
</html>"),
    DymolaStoredErrors);
end BaseWall;
