within ThermalSeparation.PressureLoss.RandomPackedColumn;
model NominalLinear "linear pressure drop with nominal values"

  extends
    ThermalSeparation.PressureLoss.RandomPackedColumn.BasicPressureLossPacked;
  parameter SI.Pressure deltaP_nom=0.05e5 "nominal pressure drop";
  parameter SI.VolumeFlowRate Vdot_nom=0.4 "nominal volume flow rate";
  final parameter Real K=(deltaP_nom/Vdot_nom)/n;
SI.ReynoldsNumber Re[n];
  SI.WeberNumber We[n];
  SI.FroudeNumber Fr[n];
    SI.Velocity w_eff_l[n] "effective liquid velocity";
      SI.Velocity w_sup_l[n] = Vdot_l/geometry.A "superficial liquid velocity";
      SI.Velocity w_sup_v[n] = Vdot/geometry.A;
      Real Ft[n];
       SI.DynamicViscosity eta_v[n] = propsVap.eta;
      SI.DynamicViscosity eta_l[n] = propsLiq.eta;
      Real sin_theta = Modelica.Math.sin(geometry.theta*2*Modelica.Constants.pi/360);
equation
   for j in 1:n loop
      w_eff_l[j]=w_sup_l[j]/(geometry.eps*eps_liq[j]*Modelica.Math.sin(geometry.theta/180*Modelica.Constants.pi));
   Re[j] = w_eff_l[j]*geometry.S*rho_l[j]/eta_l[j];
     We[j] = max(0,w_eff_l[j]*w_eff_l[j]*rho_l[j]*geometry.S/sigma_l[j]);
   Fr[j] =max(0, w_eff_l[j]*w_eff_l[j]/(geometry.S*Modelica.Constants.g_n));
   Ft[j] = 29.12 * (We[j]*Fr[j])^0.15*geometry.S^0.359/(Re[j]^0.2*geometry.eps^0.6*(1-0.93*0.62)*sin_theta^0.3);
 end for;
if homotopyMethod.bool_dp and homotopyMethod.useHomotopy then
  for j in 1:n - 1 loop
     max(0,sign(p[j] - p[j + 1])*max(0, (abs(p[j] - p[j + 1]))))=homotopy(actual=K*Vdot[j],simplified=homotopyMethod.dp);
  end for;

  /*** outlet: half of the n-th discreet element ***/
  max(0,sign(p[n] - p[n + 1])*max(0, (abs(p[n] - p[n + 1]))))=homotopy(actual=(K/2)*Vdot[n],simplified=homotopyMethod.dp/2);

  /*** entry: half of the first discreet element ***/
  p_v_in-p[1] = homotopy(actual=Vdot_in*K/2,simplified=homotopyMethod.dp/2);

else
  for j in 1:n - 1 loop
     max(0,sign(p[j] - p[j + 1])*max(0, (abs(p[j] - p[j + 1]))))=K*Vdot[j];
  end for;

  /*** outlet: half of the n-th discreet element ***/
  max(0,sign(p[n] - p[n + 1])*max(0, (abs(p[n] - p[n + 1]))))=(K/2)*Vdot[n];

  /*** entry: half of the first discreet element ***/
  p_v_in-p[1] = Vdot_in*K/2;
end if;

  annotation (Documentation(info="<html>
<p>Pressure loss model of the form:</p>
<pre><font style=\"color: #006400; \">delta_p&nbsp;=&nbsp;K*Vdot,&nbsp;</font>
<font style=\"color: #006400; \">where&nbsp;K&nbsp;is&nbsp;constant. K is determined from the nominal volume flow rate and the nominal pressure drop of the column.</font></pre>
</html>", revisions="<html>
<pre><font style=\"color: #006400; \">Documentation last revised: 18.7.2011</font></pre>
</html>"));
end NominalLinear;
