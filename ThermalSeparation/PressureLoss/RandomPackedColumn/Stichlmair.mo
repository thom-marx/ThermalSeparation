within ThermalSeparation.PressureLoss.RandomPackedColumn;
model Stichlmair "Stichlmair, for random and structured packings"

  extends ThermalSeparation.PressureLoss.RandomPackedColumn.BasicPressureLossPacked;

/*** fluid particle diameters ***/

  SI.Density rho_v_av[n-1];
  Real Re_av[n-1](start=fill(40,n-1));
  Real eps_liq_av[n-1];
  Real Re_in = 38;//rho_v[1]*geometry.dp/17e-6*Vdot_in/geometry.A;
  Real Re_out= 38;//rho_v[n]*geometry.dp/17e-6*Vdot[n]/geometry.A;

  parameter SI.Height h = geometry.H/n "height of one discreet element";

protected
  Real f0[n-1] "dimensionless factor for the discrete elements";
  Real f0_in= geometry.C1/Re_in + geometry.C2 /sqrt(max(1e-7,Re_in)) + geometry.C3;
    //"dimensionless factor for the inlet";
  Real f0_out = geometry.C1/Re_out + geometry.C2 /sqrt(max(1e-7,Re_out)) + geometry.C3;
   // "dimensionless factor for the outlet";
  Real f0_dash[n-1] "dimensionless factor for the discrete elements";
  Real f0_dash_in= f0_in * ((1-geometry.eps*(1-eps_liq[1]/geometry.eps))/(1-geometry.eps))^(c_in/3);
  Real f0_dash_out(start=1)= f0_out * ((1-geometry.eps*(1-eps_liq[n]/geometry.eps))/(1-geometry.eps))^(c_out/3);
  Real eps_dash[n-1];
  Real eps_dash_in= geometry.eps-eps_liq[1];
  Real eps_dash_out= geometry.eps-eps_liq[n];
  Real c[n-1];
  Real c_in = (-geometry.C1 / Re_in - geometry.C2/(2*sqrt(max(1e-7,Re_in))))/f0_in;
  Real c_out = (-geometry.C1 / Re_out - geometry.C2/(2*sqrt(max(1e-7,Re_out))))/f0_out;
  Real dp_dash[n-1];
  Real dp_dash_in=geometry.dp*((1-geometry.eps*(1-eps_liq[1]/geometry.eps))/(1-geometry.eps))^(1/3);
  Real dp_dash_out=geometry.dp*((1-geometry.eps*(1-eps_liq[n]/geometry.eps))/(1-geometry.eps))^(1/3);

equation
for j in 1:n-1 loop
  dp_dash[j] = geometry.dp*((1-geometry.eps*(1-eps_liq_av[j]/geometry.eps))/(1-geometry.eps))^(1/3);
  eps_dash[j] = geometry.eps-eps_liq_av[j];
  f0[j] = geometry.C1/Re_av[j] + geometry.C2 /sqrt(max(1e-7,Re_av[j])) + geometry.C3;
  f0_dash[j] = f0[j] * ((1-geometry.eps*(1-eps_liq_av[j]/geometry.eps))/(1-geometry.eps))^(c[j]/3);
  c[j] = (-geometry.C1 / Re_av[j] - geometry.C2/(2*sqrt(max(1e-7,Re_av[j]))))/f0[j];

    eps_liq_av[j] = (eps_liq[j]+eps_liq[j+1])/2;
  Re_av[j] =max(1, rho_v_av[j]*geometry.dp/17e-6*Vdot[j]/geometry.A);
  rho_v_av[j] = (rho_v_calc[j]+rho_v_calc[j+1])/2;
end for;

if homotopyMethod.bool_dp and homotopyMethod.useHomotopy then
  for j in 1:n-1 loop
  max(0,sign(p[j]-p[j+1])*abs(p[j]-p[j+1]))=homotopy(actual=Vdot[j]^2 /4/eps_dash[j]^4.65/dp_dash[j] *(3*f0_dash[j]*(1-eps_dash[j])*rho_v_av[j]*h)/geometry.A^2,simplified=homotopyMethod.dp);
  end for;

  /*** outlet: half of the n-th discreet element ***/

  max(0,sign(p[n]-p[n+1])*max(0,abs(p[n]-p[n+1])))=homotopy(actual=Vdot[n]^2/4/eps_dash_out^4.65/dp_dash_out*(3*f0_dash_out*(1 - eps_dash_out)*rho_v[n]*h/2)/geometry.A^2,simplified=homotopyMethod.dp/2);
  /*** entry: half of the first discreet element ***/

  p_v_in-p[1] = homotopy(actual=sign(Vdot_in)*Vdot_in*Vdot_in/geometry.A^2*3/4*f0_dash_in*(1-eps_dash_in)/eps_dash_in^4.65/rho_v[1]/dp_dash_in*h/2,simplified=homotopyMethod.dp/2);
else
  for j in 1:n-1 loop
  max(0,sign(p[j]-p[j+1])*abs(p[j]-p[j+1]))=Vdot[j]^2 /4/eps_dash[j]^4.65/dp_dash[j] *(3*f0_dash[j]*(1-eps_dash[j])*rho_v_av[j]*h)/geometry.A^2;
  end for;

  /*** outlet: half of the n-th discreet element ***/

  max(0,sign(p[n]-p[n+1])*max(0,abs(p[n]-p[n+1])))=Vdot[n]^2/4/eps_dash_out^4.65/dp_dash_out*(3*f0_dash_out*(1 - eps_dash_out)*rho_v[n]*h/2)/geometry.A^2;
  /*** entry: half of the first discreet element ***/

  p_v_in-p[1] = sign(Vdot_in)*Vdot_in*Vdot_in/geometry.A^2*3/4*f0_dash_in*(1-eps_dash_in)/eps_dash_in^4.65/rho_v[1]/dp_dash_in*h/2;
end if;
  annotation (Documentation(revisions="<html>
<pre>Documentation&nbsp;last&nbsp;revised:&nbsp;18.7.2011</pre>
</html>", info="<html>
<pre><font style=\"color: #006400; \">[1]&nbsp;Stichlmair,&nbsp;Bravo,&nbsp;Fair:&nbsp;General&nbsp;model&nbsp;for&nbsp;prediciton&nbsp;of&nbsp;pressure&nbsp;drop&nbsp;and&nbsp;capacity&nbsp;of&nbsp;countercurrent&nbsp;gas/liquid&nbsp;packed&nbsp;columns,&nbsp;Gas&nbsp;Separation&nbsp;&AMP;AMP;&nbsp;Purification,&nbsp;1989,&nbsp;Vol&nbsp;3,&nbsp;p.&nbsp;19-&nbsp;28</font></pre>
</html>"));
end Stichlmair;
