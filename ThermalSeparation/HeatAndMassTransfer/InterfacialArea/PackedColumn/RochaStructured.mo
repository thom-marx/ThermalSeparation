within ThermalSeparation.HeatAndMassTransfer.InterfacialArea.PackedColumn;
model RochaStructured "structued packing - Rocha et al."

  extends BasePacked;
  SI.ReynoldsNumber Re[n];
  SI.WeberNumber We[n];
  SI.FroudeNumber Fr[n];
  SI.Velocity w_eff[n] "effective liquid velocity";
   Real omega = 1;//0.5*tanh(0.08*(time-200))+0.5;
equation
  for j in 1:n loop
   w_eff[j]=w_sup[j]/(geometry.eps*eps_liq[j]*Modelica.Math.sin(geometry.theta/180*Modelica.Constants.pi));
   Re[j] = w_eff[j]*geometry.S*rho[j]/eta[j];
   We[j] = max(0,w_eff[j]*w_eff[j]*rho[j]*geometry.S/sigma[j]);
   Fr[j] =max(0, w_eff[j]*w_eff[j]/(geometry.S*Modelica.Constants.g_n));
  frac[j] = geometry.F_SE * 29.12*(We[j]*Fr[j])^0.15*geometry.S^0.359/(Re[j]^0.2*geometry.eps^0.6 * (1-0.93*geometry.cos_gamma)*(Modelica.Math.sin(geometry.theta/180*Modelica.Constants.pi))^0.3);
  a[j]=geometry.a*(1-omega)+omega*frickel*geometry.a*max(0.1,frac[j]);
 //a[j]=geometry.a;
  end for;
  annotation (Documentation(info="<html>
<p>The correlation from Rocha is suitable for structured packings. The equation was taken from the following publication:</p>
<p>Rocha, Bravo, Fair: Distillation Columns Containing Structured Packings: A comprehensive model for their performance. 2. Mass Transfer Model, Ind. Eng. Chem. Res., 1996, <i>35</i>, 1660-1667</p>
</html>", revisions="<html>
<p>Documentation last revised: 18.7.2011</p>
</html>"));
end RochaStructured;
