within ThermalSeparation.Media.Correlations;
package HenryCoefficient
  "this package contains different models to use Henry's law"

model BaseHenry

  parameter Integer nS(min=1);
  input SI.MoleFraction x_l[nS];
  input SI.Temperature T;
  output Real He[nS];

end BaseHenry;

  model Exponential "exponential approach for temperature dependency"
    extends BaseHenry;
    parameter Boolean henry_temp;
    parameter SI.Pressure henry_H[nS] "Henry coefficient";
    parameter Real henry_C[nS] "constant to calculate temperature dependency";
    parameter SI.Temperature henry_T_norm=298 "norm temperature";
  equation

   for i in 1:nS loop
     He[i] = if henry_temp then henry_H[i] * exp(henry_C[i]*(1/T - 1/henry_T_norm)) else henry_H[i];
   end for;

  end Exponential;
end HenryCoefficient;
