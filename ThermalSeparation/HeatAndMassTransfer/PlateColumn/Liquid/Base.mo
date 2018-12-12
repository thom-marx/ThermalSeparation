within ThermalSeparation.HeatAndMassTransfer.PlateColumn.Liquid;
partial model Base
  extends ThermalSeparation.HeatAndMassTransfer.PlateColumn.Base;
input SI.Height h[n] "height of the two-phase regime on the tray";
input Real eps_liq_2ph[n]
    "fraction of liquid in the two-phase regime on the tray";
  input SI.DynamicViscosity eta_comp[n,nS]
    "viscosity of each component separately";

replaceable package MediumLiquid = 
    Media.BaseMediumLiquid;

  MediumLiquid.DiffusionCoefficient[n] diffCoeff(T=T, p=p, x=x, eta=eta_comp);

protected
  SI.DiffusionCoefficient D[n,a]=diffCoeff.D;
  SI.Velocity w_sup[n] "superficial vapour velocity";
equation
           for j in 1:n loop
  w_sup[j] = max(1e-10,Vdot_l[j]/geometry.A);
  end for;

end Base;
