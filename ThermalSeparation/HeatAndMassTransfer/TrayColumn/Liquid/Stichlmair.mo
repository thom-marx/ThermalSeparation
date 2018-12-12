within ThermalSeparation.HeatAndMassTransfer.TrayColumn.Liquid;
model Stichlmair
 extends BaseLiqMT;

equation
  for j in 1:n loop
    for i in 1:n_k loop
    k[j,i] = max(1e-5,sqrt(4/Modelica.Constants.pi .* D[j,i] * w_sup[j]/max(1e-5,abs(h[j]))/max(1e-5,abs((1-eps_liq_2ph[j])))));
    //k_2[j,:] = sqrt(4/Modelica.Constants.pi .* D[j,:] * w_sup[j]);
    end for;
  end for;

                                                  annotation(choicesAllMatching,
      Documentation(info="<html>
<pre><font style=\"color: #006400; \">Stichlmair:&nbsp;Grundlagen&nbsp;der&nbsp;Dimensionierung&nbsp;des&nbsp;Gas/Fl&uuml;ssigkeit-Kontaktapparates&nbsp;Bodenkolonne,&nbsp;Weinheim,&nbsp;NY:&nbsp;Verlag&nbsp;Chemie,&nbsp;1978,&nbsp;p.&nbsp;142</font></pre>
</html>"));
end Stichlmair;
