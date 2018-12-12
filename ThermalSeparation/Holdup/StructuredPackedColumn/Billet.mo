within ThermalSeparation.Holdup.StructuredPackedColumn;
model Billet "correlation from Billet"
  extends ThermalSeparation.Holdup.StructuredPackedColumn.BaseHoldup;

equation
    for j in 1:n loop
    eps_liq[j] = (hu_stat[j]+ hu_dyn[j]);
      hu_stat[j] = 0.002;//0.033 * exp(-0.22 * Modelica.Constants.g_n * rho[j]/(sigma[j] * geometry.a * geometry.a));

//hu_dyn[j] = (max(1e-10,3*1e-3 * Vdot[j]/geometry.A/(rho[j]*geometry.a*Modelica.Constants.g_n *Modelica.Math.sin(geometry.theta*2*Modelica.Constants.pi/360))))^(1/3)*geometry.a;
   Vdot[j] =max(0,sign(hu_dyn[j])*(abs(hu_dyn[j])) ^3 * Modelica.Constants.g_n * geometry.A * rho[j] *Modelica.Math.sin(geometry.theta*2*Modelica.Constants.pi/360)/(3*geometry.a^2 * 1e-3));
end for;
  annotation (Documentation(info="<html>
</html>"));
end Billet;
