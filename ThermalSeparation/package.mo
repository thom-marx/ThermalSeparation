within ;
package ThermalSeparation "Thermal Separation Library"
constant Integer inter[:] = {1,3,6,10,15, 21, 28, 36, 45}
  "number of binary interaction values depending on the number of substances";
import SI = Modelica.SIunits;



annotation (preferedView="info",
  version="0.2",
uses(Modelica(version="3.2.3")),
     Icon(graphics={Ellipse(
        extent={{-70,70},{70,-70}},
        lineColor={0,0,0},
        fillPattern=FillPattern.Solid,
        fillColor={255,0,0}), Line(
        points={{-100,-70},{-70,-70},{70,70},{100,70}},
        color={0,0,0},
        smooth=Smooth.None)}),
    Documentation(info="<html>
<p>This library is intended to describe the following separation units:</p>
<ul>
<li>Rectification column</li>
<li>Absorption column </li>
</ul>
<p><br>For getting started please take a look at the<a href=\"Modelica://ThermalSeparation.UsersGuide\"> User&apos;s Guide</a>.</p>
<p><img src=\"../images/tt_tuhh_logo.gif\"/></p>
<p><br>Technische Universit&auml;t Hamburg</p>
<p>Institut f&uuml;r Technische Thermodynamik </p>
<p>Denickestra&szlig;e 17 </p>
<p>D-21073 Hamburg </p>
<p>Germany </p>
<p><a href=\"http://www.tu-harburg.de/tt\">www.tu-harburg.de/tt</a></p>
<p>Copyright &copy; 2017, Institute of Engineering ThermoDynamics (Hamburg University of Technology (TUHH)).</p>
<p><i>This Modelica package is <b>free</b> software; it can be redistributed and/or modified under the terms of the <b>Modelica license</b>, see the license conditions and the accompanying <b>disclaimer</b> <a href=\"../ModelicaLicense2.html\">here</a>.</i> </p>
</html>", revisions="<html>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"1\"><tr>
<td><p align=\"center\"><h4>created by </h4></p></td>
<td><p>Karin Dietl &AMP; Andreas Joos &AMP; Kai Wellner &AMP; Thomas Marx-Schubach </p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>creation date </h4></p></td>
<td><p>01.01.2017</p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>revised by </h4></p></td>
<td><p>nobody so far</p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>last revision </h4></p></td>
<td><p>this is an alpha version </p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>based on </h4></p></td>
<td></td>
</tr>
<tr>
<td><p align=\"center\"><br><h4>requires </h4></p></td>
<td><p>MSL 3.2.2</p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>tested on </h4></p></td>
<td><p>Dymola 2017 FD01</p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>contains sensitive data </h4></p></td>
<td><p>yes</p></td>
</tr>
</table>
</html>"),
  conversion(noneFromVersion="0.1"));
end ThermalSeparation;
