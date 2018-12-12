within ThermalSeparation.HeatAndMassTransfer.SprayColumn;
package Liquid
  model Constant
   // extends ThermalSeparation.FilmModel.BaseClasses.SprayColumn.BaseMSLiquid;
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
</html>",   revisions="<html>
<pre>Documentation&nbsp;last&nbsp;revised:&nbsp;18.7.2011</pre>
</html>"));
  end Constant;

  partial model BaseLiqMT "base model for liquid mass transfer coefficient"
    parameter Integer n(min=1) annotation(Dialog(enable=false));
    parameter Integer n_k
      "number of calculated mass transfer coefficients - nSL or k";
    output ThermalSeparation.Units.CoefficentOfMassTransfer k[n,n_k];
  equation

  end BaseLiqMT;
end Liquid;
