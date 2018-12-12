within ThermalSeparation.FilmModel.BaseClasses.TrayColumn;
model BaseMSVapour
  "base model for Maxwell-Stefan R-matrix (vapour side) in a tray column"
  import ThermalSeparation;
  extends ThermalSeparation.FilmModel.BaseClasses.TrayColumn.BaseMSTrayColumn;

input SI.Height h[n] "height of the two-phase regime on the tray";
input Real eps_liq_2ph[n]
    "fraction of liquid in the two-phase regime on the tray";

  replaceable package MediumVapour = 
    Media.BaseMediumVapour;

  MediumVapour.DiffusionCoefficient[n] diffCoeff(T=T, p=p);

   replaceable model MassTransferCoeff = 
     ThermalSeparation.HeatAndMassTransfer.TrayColumn.Vapour.BaseVapMT;
 MassTransferCoeff massTransferCoeff(n=n, n_k=a,D=D, w_sup=w_sup, eps_liq_2ph=eps_liq_2ph, h=h,
 redeclare record Geometry =  Geometry);

protected
  SI.DiffusionCoefficient D[n,a]=diffCoeff.D;

     SI.Velocity w_sup[n] "superficial vapour velocity";
equation
    for j in 1:n loop

    w_sup[j] = max(1e-10,Vdot_v[j]/geometry.A);
    end for;
      k_2=massTransferCoeff.k;
end BaseMSVapour;
