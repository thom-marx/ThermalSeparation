within ThermalSeparation.Media.Correlations;
package SaturationFugacityCoefficient


  annotation (Documentation(info="<html>
<p><h4><font color=\"#008000\">Liquid Phase</font></h4></p>
<p><u><font style=\"color: #008000; \">Saturation Fugacity Coefficient:</font></u></p>
<p><i><b>phi_sat</b></i> is given by the virial equation of state and depends on the systems temperature, the vapor pressure of the substance and the second virial coefficient for the pure component.</p>
<p>The second virial coefficient <i><b>Bii</b></i> are calculated with the Tsonopoulos relations, which demand specific constants. These Tsonopoulos constants are functions of the critical data and the dipole moment. For the extisting media models the relations are implemented, for further models it is necessary to implement a vector (<i><b>eq_Tsonopoulos</b></i>) for the substance specific equation. The Tsonopoulos relations and the equations for the constants can be found in (Lit [1], 4.15).</p>
<p><br/><u><font style=\"color: #008000; \">References:</font></u></p>
<p>[1] Poling, Prausnitz, O&apos;Connell : The Properties of Gases &AMP; Liquids 5th Edition [4.13, 5.10]</p>
</html>"));
end SaturationFugacityCoefficient;
