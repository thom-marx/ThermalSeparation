within ThermalSeparation.HeatAndMassTransfer.InterfacialArea.PackedColumn;
model Onda "correlation from Onda, for random packings"
//aus: Mersmann: Thermische Verfahrenstechnik, p. 337
  extends BasePacked;
equation
  for j in 1:n loop
    frac[j]=1-exp(-1.45*(geometry.sigma_crit/sigma[j])^0.75 * (w_sup[j] *rho[j]/(geometry.a * eta[j]))^0.1 * (w_sup[j]^2 * geometry.a/Modelica.Constants.g_n)^(-0.05) * (w_sup[j]^2 * rho[j]/(sigma[j]*geometry.a))^0.2);
  //a[j]=if time <1 then geometry.a*0.5 else geometry.a*max(0.1,frac[j]);
  a[j]= frickel*geometry.a*max(0.1,frac[j]);
  end for;
  annotation (Documentation(info="<html>
<p>This class calculates the interfacial area using the correlation from Onda. This approach is suitable for random packings. The equation itself was copied from: Mersmann, Kind, Stichlmair: Thermische Verfahrenstechnik, Springer 2005, p. 337.</p>
</html>", revisions="<html>
<p>Documentation last edited: 18.7.2011</p>
</html>"));
end Onda;
