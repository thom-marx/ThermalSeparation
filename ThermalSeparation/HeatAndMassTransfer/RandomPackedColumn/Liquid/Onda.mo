within ThermalSeparation.HeatAndMassTransfer.RandomPackedColumn.Liquid;
model Onda "correlation from Onda, for random packings"
//extends ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSLiquid;
  extends BaseLiqMT;
  parameter ThermalSeparation.Units.CoefficentOfMassTransfer k_l_0=
                                                        1e-4
    "liquid mass transfer coefficient, if velocity is 0";

  ThermalSeparation.Units.VolumetricArea a_wetted[
                                     n];

//variables used in a_wetted
  Real aA[n];
  Real aB[n];
  Real aC[n];
  Real aD[n];
//variables used in k_l
protected
  Real klA[n];
  Real klB[n];
  Real klC[n,n_k];
  Real klD[n];
  Boolean liquidExists[n](each start=false) "true if liquid is on the stage";

equation
  // Werte belegen, damit das System nicht singulär wird:
//   Re[:]=fill(0,n);
//   Sc[:,:]=fill(0,n,a);
//   Sh[:,:]=fill(0,n,a);
  for j in 1:n loop
    a_wetted[j] = max(1e1,geometry.a *(1- exp(-1.45*(aA[j])^0.75 * (aB[j])^0.1 * (aC[j])^(-0.05) * (aD[j])^0.2)));
      aA[j] = geometry.sigma_crit/sigma[j];
      aB[j] = w_sup[j]*rho[j]/geometry.a/eta[j];
      aC[j] = max(1e-10,w_sup[j]*w_sup[j] * geometry.a / Modelica.Constants.g_n);
      aD[j] = w_sup[j]*w_sup[j]*rho[j]/sigma[j]/geometry.a;
    when w_sup[j]>0 then
      liquidExists[j] = true;
    end when;
    k[j,:] = if liquidExists[j] then 0.0051*(klA[j])^(-1/3) * (klB[j])^(-2/3) * klC[j,:].^(-1/2) * (klD[j])^(1/5) else fill(k_l_0,n_k);
      klA[j] = rho[j]/(Modelica.Constants.g_n*eta[j]);
      klB[j] = max(1e-10,rho[j]*w_sup[j]/(a_wetted[j]*eta[j]));
      klC[j,:] = eta[j]./(rho[j].*D[j,:]);
      klD[j] = geometry.a * geometry.d_char;
end for;
  annotation (Documentation(info="<html>
<p>Correlation from Onda [1] is used in order to calculate the binary mass transfer coefficients.</p>
<p>[1] Onda et al.: Mass transfer coefficients between gas and liquid phases in packed columns, Journal of Chem. Eng. of Japan, Vol. 1 (1968), p. 56-62</p>
</html>", revisions="<html>
<pre>Documentation&nbsp;last&nbsp;revised:&nbsp;18.7.2011</pre>
</html>"));
end Onda;
