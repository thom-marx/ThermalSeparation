within ThermalSeparation.FilmModel.BaseClasses.StructuredPackedColumn;
model BaseMSVapour
  "base model for Maxwell-Stefan R-matrix (vapour side) in a packed column"
  extends
    ThermalSeparation.FilmModel.BaseClasses.StructuredPackedColumn.BaseMSPackedColumn;

  input SI.SurfaceTension sigma[n];

  replaceable package MediumVapour = 
    Media.BaseMediumVapour;

  MediumVapour.DiffusionCoefficient[n] diffCoeff(T=T, p=p);
   replaceable model MassTransferCoeff = 
     ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Vapour.BaseVapMT;
 MassTransferCoeff massTransferCoeff(n=n, n_k=a,D=D, eta=eta, sigma=sigma, rho=rho, w_sup_l=w_sup_l, w_sup_v=w_sup_v, eps_liq=eps_liq,
 redeclare record Geometry =  Geometry);

protected
  SI.DiffusionCoefficient D[n,a]=diffCoeff.D;
     SI.Velocity w_sup_l[n] "superficial liquid velocity";
     SI.Velocity w_sup_v[n] "superficial vapour velocity";
equation
    for j in 1:n loop
  w_sup_l[j] = max(1e-10,Vdot_l[j]/geometry.A);
    w_sup_v[j] = max(1e-10,Vdot_v[j]/geometry.A);
    end for;
    k_2=massTransferCoeff.k;
end BaseMSVapour;
