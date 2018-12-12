within ThermalSeparation.HeatAndMassTransfer.PlateColumn.Vapour;
model Base
  import ThermalSeparation;
  extends ThermalSeparation.HeatAndMassTransfer.PlateColumn.Base;

input SI.Height h[n] "height of the two-phase regime on the tray";
input Real eps_liq_2ph[n]
    "fraction of liquid in the two-phase regime on the tray";

  replaceable package MediumVapour = 
      ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2_CO2_SO2_HCl_HF 
                                                            constrainedby
    Media.BaseMediumVapour;

  MediumVapour.DiffusionCoefficient[n] diffCoeff(T=T, p=p);

protected
  SI.DiffusionCoefficient D[n,a]=diffCoeff.D;

     SI.Velocity w_sup[n] "superficial vapour velocity";
equation
    for j in 1:n loop

    w_sup[j] = max(1e-10,Vdot_v[j]/geometry.A);
  end for;
end Base;
