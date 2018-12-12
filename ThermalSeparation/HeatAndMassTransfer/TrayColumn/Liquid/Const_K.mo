within ThermalSeparation.HeatAndMassTransfer.TrayColumn.Liquid;
model Const_K
extends BaseLiqMT;
  parameter ThermalSeparation.Units.CoefficentOfMassTransfer k_l_const=
                                                          1e-3
    "mass transfer coefficient for liquid";

//Wesselingh: Mass Transfer in Multicomponent Mixtures, p. 52
//gases: 1e-1 to 1e-2 (for gases in pores, the value can be smaller, for example 5e-3)
//liquids: 1e-4 to 1e-5 (for liquids in pores, the value can be smaller, for example 1e-6)
  //1e-4, 1e-2
equation

  for j in 1:n loop
    k[j,:] = fill(k_l_const,n_k);
  end for;
  annotation (Documentation(info="<html>
<p>Constant value for liquid mass transfer coefficient.</p>
</html>"));
end Const_K;
