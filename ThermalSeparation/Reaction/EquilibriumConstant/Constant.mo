within ThermalSeparation.Reaction.EquilibriumConstant;
model Constant "constant value for K_eq"
  extends ThermalSeparation.Reaction.EquilibriumConstant.BaseK;
  parameter Real[nR] K_const = fill(100,nR);
  //0.0139
equation

    K[:] = K_const;

end Constant;
