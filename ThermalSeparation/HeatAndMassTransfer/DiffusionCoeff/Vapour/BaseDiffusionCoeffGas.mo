within ThermalSeparation.HeatAndMassTransfer.DiffusionCoeff.Vapour;
partial model BaseDiffusionCoeffGas
parameter Integer n=1;
parameter Integer nS=2;
parameter Integer a = 1;

input SI.Temperature T[n];
input SI.Pressure p[n];

replaceable package MediumVapour = 
    Media.BaseMediumVapour;
  //Diff.koeff. z.B.: D12, D13, D14, D23, D24, D34
  output SI.DiffusionCoefficient D[n,a];

end BaseDiffusionCoeffGas;
