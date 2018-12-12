within ThermalSeparation.HeatAndMassTransfer.PlateColumn.Vapour;
model Const_K
  extends Base;

  parameter ThermalSeparation.Units.CoefficentOfMassTransfer k_v_const=
                                                          1e-3
    "mass transfer coefficient for gas";

//Wesselingh: Mass Transfer in Multicomponent Mixtures, p. 52
//gases: 1e-1 to 1e-2 (for gases in pores, the value can be smaller, for example 5e-3)
//liquids: 1e-4 to 1e-5 (for liquids in pores, the value can be smaller, for example 1e-6)
  //1e-4, 1e-2
equation
    // Werte belegen, damit das System nicht singulär wird:
//   Re[:]=fill(0,n);
//   Sc[:,:]=fill(0,n,a);
//   Sh[:,:]=fill(0,n,a);
  for j in 1:n loop

    k_2[j,:] = fill(k_v_const,a);
  end for;

end Const_K;
