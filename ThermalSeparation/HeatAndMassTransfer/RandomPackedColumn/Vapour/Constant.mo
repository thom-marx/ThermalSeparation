within ThermalSeparation.HeatAndMassTransfer.RandomPackedColumn.Vapour;
model Constant "constant value for k"
 // extends ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSVapour;
extends BaseVapMT;
  parameter ThermalSeparation.Units.CoefficentOfMassTransfer k_v_const=
                                                          1e-2
    "mass transfer coefficient for gas";

//Wesselingh: Mass Transfer in Multicomponent Mixtures, p. 52
//gases: 1e-1 to 1e-2 (for gases in pores, the value can be smaller, for example 5e-3)
//liquids: 1e-4 to 1e-5 (for liquids in pores, the value can be smaller, for example 1e-6)
equation

    // Werte belegen, damit das System nicht singulär wird:
//   Re[:]=fill(0,n);
//   Sc[:,:]=fill(0,n,a);
//   Sh[:,:]=fill(0,n,a);
    for j in 1:n loop
    k[j,:] = fill(k_v_const,n_k);
  end for;
  annotation (Documentation(info="<html>
<p>Constant value for binary mass transfer coefficients.</p>
</html>", revisions="<html>
<pre>Documentation&nbsp;last&nbsp;revised:&nbsp;18.7.2011</pre>
</html>"));
end Constant;
