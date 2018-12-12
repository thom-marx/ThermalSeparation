within ThermalSeparation.Media.Correlations.DiffusionCoefficient.Vapour;
model Const_Gas "constant diffusion coefficient for gas"
  extends
    ThermalSeparation.Media.Correlations.DiffusionCoefficient.Vapour.BaseDiffusionCoeffGas;
  parameter SI.DiffusionCoefficient D_const_gas = 1e-5;
equation
  D=fill(D_const_gas,a);
  D_matrix = fill(D_const_gas,nS,nS);
  annotation (Documentation(info="<html>
<p>The binary diffusion coefficient is set constant. The default value is 1e-5 which is a good starting point.</p>
</html>"));
end Const_Gas;
