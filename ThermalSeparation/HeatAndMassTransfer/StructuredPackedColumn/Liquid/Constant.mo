within ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Liquid;
model Constant
 // extends ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSLiquid;
extends BaseLiqMT;
  parameter ThermalSeparation.Units.CoefficentOfMassTransfer k_l_const=
                                                          1e-4
    "mass transfer coefficient for liquid";

//Wesselingh: Mass Transfer in Multicomponent Mixtures, p. 52
//gases: 1e-1 to 1e-2 (for gases in pores, the value can be smaller, for example 5e-3)
//liquids: 1e-4 to 1e-5 (for liquids in pores, the value can be smaller, for example 1e-6)
equation

    for j in 1:n loop
    k[j,:] = fill(k_l_const,n_k);
  end for;
  annotation (Documentation(info="<html>
<p>Constant value for binary mass transfer coefficients.</p>
</html>", revisions="<html>
<pre>Documentation&nbsp;last&nbsp;revised:&nbsp;18.7.2011</pre>
</html>"));
end Constant;
