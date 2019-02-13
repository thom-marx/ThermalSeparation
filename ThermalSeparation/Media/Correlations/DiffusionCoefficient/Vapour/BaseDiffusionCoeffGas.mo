within ThermalSeparation.Media.Correlations.DiffusionCoefficient.Vapour;
partial model BaseDiffusionCoeffGas
parameter Integer nS=2;

input SI.Temperature T;
input SI.Pressure p;

  parameter SI.MolarMass MMX[nS];
  parameter SI.Volume V[nS] "diffusion volumes";

  //Diff.koeff. z.B.: D12, D13, D14, D23, D24, D34
  output SI.DiffusionCoefficient D[a](start=fill(1e-5,a));
  output SI.DiffusionCoefficient D_matrix[nS,nS];

protected
  parameter Integer aux[:] = {1,3,6,10,15, 21, 28, 36, 45};
  parameter Integer a = aux[nS-1]
    "number of binary diffusion coefficients depending on the number of substances";

end BaseDiffusionCoeffGas;
