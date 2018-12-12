within ThermalSeparation.Reaction.EquilibriumConstant;
model Simple "simple temperature dependency introduced"
  extends ThermalSeparation.Reaction.EquilibriumConstant.BaseK;
  parameter Real[nR] K0 = fill(0.5,nR)
    "equilibrium constant for temperature T0";
  parameter SI.Temperature[nR] T0 = fill(50+273,nR)
    "reference temperature for equilibrium constant";
protected
  Real A[nR];
equation
  for m in 1:nR loop
    A[m] = -Modelica.Math.log(K0[m]) * T0[m];
  end for;

    for m in 1:nR loop
      K[m] = exp(-1/T);
    end for;

end Simple;
