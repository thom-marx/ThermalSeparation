within ThermalSeparation.Holdup.StructuredPackedColumn;
model StichlmairStat "correlation from Stichlmair, for random packings"
  extends ThermalSeparation.Holdup.StructuredPackedColumn.BaseHoldup;
equation
    for j in 1:n loop

    eps_liq[j] = (hu_stat[j]+ hu_dyn[j]);
      hu_stat[j] = 0.033 * exp(-0.22 * Modelica.Constants.g_n * rho[j]/(max(1e-7,sigma[j] * geometry.a * geometry.a)));

   Vdot[j] =max(0,sign(hu_dyn[j])*((abs(hu_dyn[j]) / 0.555)^3 * Modelica.Constants.g_n * geometry.eps^4.65/ geometry.a)^0.5* geometry.A * geometry.eps *eps_liq[j]*50);
end for;
  annotation (Documentation(info="<html>
<p>The correlation for the static holdup was taken from [1] (equation 2), the correlation for the dynamic holdup was taken from [2] (equation 5).</p>
<p><br/>References:</p>
<p>[1] Engel, V.: Fluiddynamik in F&uuml;llk&ouml;rper- und Packungskolonnen f&uuml;r Gas/Fl&uuml;ssig-Systeme, Chemie Ingenieur Technik, Vol. 72, pp. 700-703</p>
<p>[2] Wagner, I et al.: Mass Transfer in Beds of Modern, High-Efficiency Random Packing, Ind. Eng. Chem. Res, Vol 36 (1997), pp. 227-237</p>
</html>"));
end StichlmairStat;
