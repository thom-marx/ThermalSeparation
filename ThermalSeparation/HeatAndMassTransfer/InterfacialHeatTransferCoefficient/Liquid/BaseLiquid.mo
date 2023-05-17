within ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Liquid;
partial model BaseLiquid
  extends ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.BaseInterfacialHTCoefficient;
    replaceable package Medium =
    Media.BaseMediumLiquid;
    input Medium.ThermodynamicProperties props[n];
equation

end BaseLiquid;
