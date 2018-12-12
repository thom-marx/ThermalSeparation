within ThermalSeparation.PressureLoss.PackedColumn;
model NominalQuadratic "quadratic pressure loss with nominal values"
//delta_p = K*Vdot^2, where K is constant
  extends ThermalSeparation.PressureLoss.PackedColumn.BasicPressureLossPacked;
  parameter SI.Pressure deltaP_nom= 0.05e5 "nominal pressure drop";
  parameter SI.VolumeFlowRate Vdot_nom = 0.1 "nominal volume flow rate";
protected
  Real K = deltaP_nom/Vdot_nom^2;
equation
for j in 1:n-1 loop
Vdot[j] = max(0,sign(p[j]-p[j+1])*sqrt(abs(p[j]-p[j+1])/(K/n)));
end for;

/*** outlet: half of the n-th discreet element ***/
Vdot[n] = max(0,sign(p[n]-p[n+1])*sqrt(max(1e-5,abs(p[n]-p[n+1])/(K/n)*2)));

/*** entry: half of the first discreet element ***/
p_v_in = p[1]+sign(Vdot_in)*Vdot_in*Vdot_in*(K/n)/2;

  annotation (Documentation(info="<html>
<p>Pressure loss model of the form:</p>
<pre><font style=\"color: #006400; \">delta_p&nbsp;=&nbsp;K*Vdot^2,&nbsp;</font>
<font style=\"color: #006400; \">where&nbsp;K&nbsp;is&nbsp;constant. K is determined from the nominal volume flow rate and the nominal pressure drop of the column.</font></pre>
</html>", revisions="<html>
<pre><font style=\"color: #006400; \">Documentation last revised: 18.7.2011</font></pre>
</html>"));
end NominalQuadratic;
