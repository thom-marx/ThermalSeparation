within ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Vapour;
model Constant "constant value"
  extends BaseVapour;
  parameter SI.CoefficientOfHeatTransfer alpha_const= 1e5;
equation
  alpha=fill(alpha_const,n);
end Constant;
