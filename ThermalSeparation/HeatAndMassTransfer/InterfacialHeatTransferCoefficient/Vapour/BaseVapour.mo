within ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Vapour;
partial model BaseVapour
  extends
    ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.BaseInterfacialHTCoefficient;
  replaceable package Medium =
    Media.BaseMediumVapour;

  input Medium.ThermodynamicProperties props[n];

end BaseVapour;
