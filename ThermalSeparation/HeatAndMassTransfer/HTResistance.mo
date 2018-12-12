within ThermalSeparation.HeatAndMassTransfer;
package HTResistance

  partial model BaseHTResistance "heat transfer between bulk and wall"

  parameter Integer n(min=1) annotation(dialog(enable=false));

  input Modelica.SIunits.Temperature T[n];
  input Modelica.SIunits.Area A[n];
  input Modelica.SIunits.HeatFlowRate Qdot[n](start=fill(1e8,n));
  input Modelica.SIunits.Pressure p[n];

  Modelica.SIunits.Temperature Twall[n];
  Real alpha[n] "heat transfer coefficient";

  end BaseHTResistance;

  model BulkBoiling "for bulk boiling in pressure vessels (natural convection)"
    extends BaseHTResistance;

  equation
    for i in 1:n loop
      Qdot[i] = alpha[i] * A[i] * (Twall[i]-T[i]);

      alpha[i]=1.95*(abs(Qdot[i])/A[i])^0.72*(p[i]/1e5)^0.24; // Stephan, Baehr
    end for;
  end BulkBoiling;

  model ConstantAlpha "heat transfer resistance with const. HT coefficient"
    extends BaseHTResistance;

  parameter Real alpha_user = 2e4 "heat transfer coefficient";

  equation
    for i in 1:n loop
      Qdot[i] = alpha[i] * A[i] * (Twall[i]-T[i]);

      alpha[i]=alpha_user;
    end for;
  end ConstantAlpha;

  model NoHTResistance "no heat transfer resistance between bulk and wall"
    extends BaseHTResistance;

  equation
    for i in 1:n loop
      alpha[i] = 1e8; // actually alpha goes to infinity

      Twall[i] = T[i];
    end for;
  end NoHTResistance;
end HTResistance;
