within ThermalSeparation.FilmModel.BaseClasses.SprayColumn;
model BaseMSLiquid
  "base model for Maxwell-Stefan R-matrix (liquid side) in a spray column"
  extends ThermalSeparation.FilmModel.BaseClasses.SprayColumn.BaseMSSprayColumn;
     replaceable model MassTransferCoeff = 
     ThermalSeparation.HeatAndMassTransfer.SprayColumn.Liquid.BaseLiqMT;
 MassTransferCoeff massTransferCoeff(n=n, n_k=a);

equation
  k_2=massTransferCoeff.k;
end BaseMSLiquid;
