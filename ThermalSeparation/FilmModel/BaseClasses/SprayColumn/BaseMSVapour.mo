within ThermalSeparation.FilmModel.BaseClasses.SprayColumn;
model BaseMSVapour
  "base model for Maxwell-Stefan R-matrix (vapour side) in a spray column"
  extends ThermalSeparation.FilmModel.BaseClasses.SprayColumn.BaseMSSprayColumn;
     replaceable model MassTransferCoeff = 
     ThermalSeparation.HeatAndMassTransfer.SprayColumn.Vapour.BaseVapMT;
 MassTransferCoeff massTransferCoeff(n=n, n_k=a);

equation
  k_2=massTransferCoeff.k;
end BaseMSVapour;
