within ThermalSeparation.Holdup.RandomPackedColumn;
model NoOutflow "sump: no liquid flowing out"
  extends ThermalSeparation.Holdup.RandomPackedColumn.BaseHoldup;
  //Rocha, Bravo, Fair: Distillation Columns Containing Structured Packings: A Comprehensive Model for Their Performance 1. Hydraulic Models, Ind. Eng. Chem. Res., 1993, 32, 641-651
  //only for structured packings
parameter Boolean hu_stat_const = true
    "true, if a constant value for static holdup is to be used";
parameter Real hu_stat_user = 0.02 annotation(Dialog(enable=hu_stat_const));
parameter Real scf_stat = 0.01
    "surface correction factor, depending on package surface";
parameter Real frickel_dyn = 1 "factor to correct the dynamic holdup";
parameter Boolean use_rho_l_nom=true;
parameter Boolean use_rho_v_nom=true;
parameter Boolean use_sigma_nom=true;
parameter Boolean use_eta_nom=true;
parameter SI.Density rho_l_nom = 1100;
SI.Density rho_l[n] = if use_rho_l_nom then fill(rho_l_nom,n) else rho;
parameter SI.Density rho_v_nom = 1.2;
SI.Density rho_v2[n] = if use_rho_v_nom then fill(rho_v_nom,n) else rho_v;
parameter Real sigma_nom = 0.072;
Real sigma_l[n]= if use_sigma_nom then fill(sigma_nom,n) else sigma;
parameter SI.DynamicViscosity eta_nom = 1e-3;
SI.DynamicViscosity eta_l[n] = if use_eta_nom then fill(eta_nom,n) else eta;

protected
  Real a;
  Real b[n];
  Real sin_theta = Modelica.Math.sin(geometry.theta*2*Modelica.Constants.pi/360);
  Real c[n];
  parameter SI.Length deltaZ = geometry.H/n;
  Real K1[n] = Modelica.Constants.g_n * (rho_l-rho_v)/1025
    "dimensionless constant from maximum capacity data";

equation
  a= 4/geometry.S *7.2453*0.76*5.64*geometry.S^0.397/(geometry.eps^0.6*(1-0.93*geometry.cos_gamma)*(sin_theta)^0.3);

    for j in 1:n loop
    eps_liq[j] = (hu_stat[j]+ hu_dyn[j]);

    if hu_stat_const then
      hu_stat[j] = hu_stat_user;
    else
      hu_stat[j] = (2*sigma_l[j]*(1-geometry.cos_gamma)/(rho_l[j]*Modelica.Constants.g_n*(1-rho_v2[j]/rho_l[j])*sin_theta))^0.5 *4/geometry.S*0.01;
    end if;

//der Term mit K1 wird vernachlässigt, der Term wird sowieso erst wichtig, wenn deltaP/deltaZ nahe am Wert von (deltaP/deltaZ)_flood ist
b[j] = 1-rho_v[j]/rho_l[j];//- K1[j]*deltaP[j]/deltaZ/(rho_l[j]* Modelica.Constants.g_n);
c[j] = (rho_l[j]/(sigma_l[j]*Modelica.Constants.g_n))^0.15 * (eta_l[j]/(geometry.S*rho_l[j]))^0.2 * (3*eta_l[j])^0.5;

Vdot[j] =0;// max(0,sign(hu_dyn[j]) * (abs(hu_dyn[j])/(a*c[j]*frickel_dyn)*(rho_l[j]*geometry.eps*eps_liq[j]*Modelica.Constants.g_n*b[j]*sin_theta)^0.5)^(1/0.9) * geometry.A);

end for;
  annotation (Documentation(info="<html>
<p>The correlation was taken from [1].</p>
<p><br/>References:</p>
<p>[1]  <code><font style=\"color: #006400; \">Rocha,&nbsp;Bravo,&nbsp;Fair:&nbsp;Distillation&nbsp;Columns&nbsp;Containing&nbsp;Structured&nbsp;Packings:&nbsp;A&nbsp;Comprehensive&nbsp;Model&nbsp;for&nbsp;Their&nbsp;Performance&nbsp;1.&nbsp;Hydraulic&nbsp;Models,&nbsp;Ind.&nbsp;Eng.&nbsp;Chem.&nbsp;Res.,&nbsp;1993,&nbsp;32,&nbsp;641-651</font></code></p>
</html>"));
end NoOutflow;
