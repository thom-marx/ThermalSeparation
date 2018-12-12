within ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Liquid;
model Constant "constant value"
  extends BaseLiquid;
   parameter SI.CoefficientOfHeatTransfer alpha_const= 1e5;
   //base model for mass transfer in a packed column
equation
  alpha=fill(alpha_const,n);
end Constant;
