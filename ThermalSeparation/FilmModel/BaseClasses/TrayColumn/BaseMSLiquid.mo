within ThermalSeparation.FilmModel.BaseClasses.TrayColumn;
model BaseMSLiquid
  "base model for Maxwell-Stefan R-matrix (liquid side) in a tray column"
  extends ThermalSeparation.FilmModel.BaseClasses.TrayColumn.BaseMSTrayColumn;
input SI.Height h[n] "height of the two-phase regime on the tray";
input Real eps_liq_2ph[n]
    "fraction of liquid in the two-phase regime on the tray";
  input SI.DynamicViscosity eta_comp[n,nS]
    "viscosity of each component separately";

replaceable package MediumLiquid = 
    Media.BaseMediumLiquid;

  MediumLiquid.DiffusionCoefficient[n] diffCoeff(T=T, p=p, x=x, eta=eta_comp);

   replaceable model MassTransferCoeff = 
     ThermalSeparation.HeatAndMassTransfer.TrayColumn.Liquid.BaseLiqMT;
 MassTransferCoeff massTransferCoeff(n=n, n_k=a,D=D, w_sup=w_sup, eps_liq_2ph=eps_liq_2ph, h=h,
 redeclare record Geometry =  Geometry);

protected
  SI.DiffusionCoefficient D[n,a]=diffCoeff.D;
  SI.Velocity w_sup[n] "superficial vapour velocity";
equation
           for j in 1:n loop
  w_sup[j] = max(1e-10,Vdot_l[j]/geometry.A);
  end for;
  k_2=massTransferCoeff.k;
end BaseMSLiquid;
