within ThermalSeparation.HeatAndMassTransfer.DiffusionCoeff.Vapour;
model Const_Gas "constant diffusion coefficient for gas"
  extends
    ThermalSeparation.HeatAndMassTransfer.DiffusionCoeff.Vapour.BaseDiffusionCoeffGas;
  parameter SI.DiffusionCoefficient D_const_gas = 1e-5;
equation
  D=fill(D_const_gas,n,a);
  annotation (Documentation(info="<html>
<p>The binary diffusion coefficient is set constant. The default value is 1e-5 which is a good starting point.</p>
</html>"));
end Const_Gas;
