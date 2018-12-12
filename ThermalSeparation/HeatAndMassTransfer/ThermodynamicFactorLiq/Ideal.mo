within ThermalSeparation.HeatAndMassTransfer.ThermodynamicFactorLiq;
model Ideal
  extends BaseThermodynamicFactor;
equation
  for j in 1:n loop
    Gamma[j,:,:]=diagonal(ones(nS));
  end for;

end Ideal;
