within ThermalSeparation.HeatAndMassTransfer.InterfacialArea;
partial model Base
  parameter Integer n(min=1);
  input SI.VolumeFlowRate Vdot_l[n];
  input SI.DynamicViscosity eta[n] "liquid viscosity";
  input SI.Density rho[n] "liquid density";
  input SI.SurfaceTension sigma[n] "liquid surface tension";
  input Real eps_liq[n];
  output Units.VolumetricArea a[n] "specific heat and mass transfer area";

end Base;
