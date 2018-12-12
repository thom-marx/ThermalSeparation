within ThermalSeparation.FilmModel.BaseClasses.StructuredPackedColumn;
model BaseMSLiquid
  "base model for Maxwell-Stefan R-matrix (liquid side) in a packed column"

 extends
    ThermalSeparation.FilmModel.BaseClasses.StructuredPackedColumn.BaseMSPackedColumn;

 input SI.SurfaceTension sigma[n];
 input SI.DynamicViscosity eta_comp[n,nS]
    "viscosity of each component separately";
 input Modelica.SIunits.Concentration c_star[n,nS]
    "concentration at phase boundary";

   replaceable package MediumLiquid = 
   Media.BaseMediumLiquid;

  MediumLiquid.DiffusionCoefficient[n] diffCoeff(
    T=T,
    p=p,
    x=x,
    eta=eta_comp);

  SI.Velocity w_sup[n] "superficial liquid velocity";
 replaceable model MassTransferCoeff = 
     ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Liquid.BaseLiqMT;
 MassTransferCoeff massTransferCoeff(n=n, n_k=a,D=D, eta=eta, sigma=sigma, rho=rho, w_sup=w_sup, eps_liq=eps_liq,
 redeclare record Geometry =  Geometry);

protected
 SI.DiffusionCoefficient D[n,a]=diffCoeff.D;
equation
  for j in 1:n loop
 w_sup[j] = max(1e-10,Vdot_l[j]/geometry.A);
  end for;
  k_2=massTransferCoeff.k;

end BaseMSLiquid;
