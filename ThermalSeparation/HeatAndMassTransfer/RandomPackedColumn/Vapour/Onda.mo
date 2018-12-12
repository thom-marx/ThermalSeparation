within ThermalSeparation.HeatAndMassTransfer.RandomPackedColumn.Vapour;
model Onda "correlation from Onda, for random packings"
 // extends ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSVapour;
extends BaseVapMT;
    parameter ThermalSeparation.Units.CoefficentOfMassTransfer k_l_0=
                                                          1e-4
    "liquid mass transfer coefficient, if velocity is 0";

equation
  for j in 1:n loop
    k[j,:] =geometry.C * (rho[j]*w_sup_v[j]/geometry.a/eta[j])^0.7 * (eta[j]/sigma[j]./D[j,:]).^(1/3) * (geometry.a*geometry.d_char)^(-2) .* (geometry.a .* D[j,:]);
  end for;

  annotation (Documentation(info="<html>
<p>Correlation from Onda [1] is used in order to calculate the binary mass transfer coefficients</p>
<p>[1] Onda et al.: Mass transfer coefficients between gas and liquid phases in packed columns, Journal of Chem. Eng. of Japan, Vol. 1 (1968), p. 56-62</p>
</html>", revisions="<html>
<pre>Documentation&nbsp;last&nbsp;revised:&nbsp;18.7.2011</pre>
</html>"));
end Onda;
