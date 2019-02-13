within ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient;
partial model BaseInterfacialHTCoefficient
  "base class for vapour and liquid interfacial heat transfer coefficient"

  parameter Integer n(min=1);
  parameter Integer nS(min=2);
  input ThermalSeparation.Units.CoefficentOfMassTransfer k_av[n]
    "averaged mass transfer coefficient";
  input SI.DiffusionCoefficient D_av[n] "averaged diffusion coefficient";
  input SI.Density rho[n];
  input SI.Concentration c[n,nS];

  output SI.CoefficientOfHeatTransfer alpha[n];
equation

end BaseInterfacialHTCoefficient;
