within ThermalSeparation.PressureLoss.StructuredPackedColumn;
model NominalQuadratic "quadratic pressure loss with nominal values"
//delta_p = K*Vdot^2, where K is constant
  extends ThermalSeparation.PressureLoss.StructuredPackedColumn.BasicPressureLossPacked;
  parameter SI.Pressure deltaP_nom= 0.05e5 "nominal pressure drop";
  parameter SI.VolumeFlowRate Vdot_nom = 0.1 "nominal volume flow rate";
protected
  Real K = deltaP_nom/Vdot_nom^2;
equation

if homotopyMethod.bool_dp and homotopyMethod.useHomotopy then
  for j in 1:n - 1 loop
     max(0,sign(p[j] - p[j + 1])*max(0, (abs(p[j] - p[j + 1]))))=homotopy(actual=K*Vdot[j]^2,simplified=homotopyMethod.dp);
  end for;

  /*** outlet: half of the n-th discreet element ***/
  max(0,sign(p[n] - p[n + 1])*max(0, (abs(p[n] - p[n + 1]))))=homotopy(actual=(K/2)*Vdot[n]^2,simplified=homotopyMethod.dp/2);

  /*** entry: half of the first discreet element ***/
  p_v_in-p[1] = homotopy(actual=Vdot_in^2*K/2,simplified=homotopyMethod.dp/2);

else
  for j in 1:n - 1 loop
     max(0,sign(p[j] - p[j + 1])*max(0, (abs(p[j] - p[j + 1]))))=K*Vdot[j]^2;
  end for;

  /*** outlet: half of the n-th discreet element ***/
  max(0,sign(p[n] - p[n + 1])*max(0, (abs(p[n] - p[n + 1]))))=(K/2)*Vdot[n]^2;

  /*** entry: half of the first discreet element ***/
  p_v_in-p[1] = Vdot_in^2*K/2;
end if;

  annotation (Documentation(info="<html>
<p>Pressure loss model of the form:</p>
<pre><font style=\"color: #006400; \">delta_p&nbsp;=&nbsp;K*Vdot^2,&nbsp;</font>
<font style=\"color: #006400; \">where&nbsp;K&nbsp;is&nbsp;constant. K is determined from the nominal volume flow rate and the nominal pressure drop of the column.</font></pre>
</html>", revisions="<html>
<pre><font style=\"color: #006400; \">Documentation last revised: 18.7.2011</font></pre>
</html>"));
end NominalQuadratic;
