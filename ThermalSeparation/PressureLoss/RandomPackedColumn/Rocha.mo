within ThermalSeparation.PressureLoss.RandomPackedColumn;
model Rocha "correlation from Rocha et al."

  extends
    ThermalSeparation.PressureLoss.RandomPackedColumn.BasicPressureLossPacked;
Real sin_theta = Modelica.Math.sin(geometry.theta*2*Modelica.Constants.pi/360);
Real A[n] = 0.177*rho_v_calc/(geometry.S * geometry.eps^2 * sin_theta^2);
Real B[n] = 88.774 * eta_v_calc /(geometry.S^2 * geometry.eps * sin_theta);
  SI.ReynoldsNumber Re[n];
  SI.WeberNumber We[n];
  SI.FroudeNumber Fr[n];
    SI.Velocity w_eff_l[n] "effective liquid velocity";
      SI.Velocity w_sup_l[n] = Vdot_l/geometry.A "superficial liquid velocity";
      SI.Velocity w_sup_v[n] = Vdot/geometry.A;
      Real Ft[n];
      parameter Real deltaP_flood = 1000;

            SI.Velocity w_sup_v_in = Vdot_in/geometry.A;

      Real C[n];

  parameter Boolean eta_const = false
    "eta is set constant (eta = eta_nom) for pressure loss calculation" annotation(Dialog(group="Use of nominal values"));
  parameter SI.Density eta_l_nom = 1e-3
                                 annotation(Dialog(enable=eta_const,group="Use of nominal values"));
  parameter SI.Density eta_v_nom = 1e-5
                                 annotation(Dialog(enable=eta_const,group="Use of nominal values"));

protected
      SI.DynamicViscosity eta_v[n] = propsVap.eta;
      SI.DynamicViscosity eta_l[n] = propsLiq.eta;
  SI.Density eta_l_calc[n] = if eta_const then fill(eta_l_nom,n) else eta_l;
  SI.Density eta_v_calc[n] = if eta_const then fill(eta_v_nom,n) else eta_v;
equation

 for j in 1:n loop
      w_eff_l[j]=w_sup_l[j]/(geometry.eps*eps_liq[j]*Modelica.Math.sin(geometry.theta/180*Modelica.Constants.pi));
   Re[j] = w_eff_l[j]*geometry.S*rho_l_calc[j]/eta_l_calc[j];
     We[j] = max(0,w_eff_l[j]*w_eff_l[j]*rho_l_calc[j]*geometry.S/sigma_l[j]);
   Fr[j] =max(0, w_eff_l[j]*w_eff_l[j]/(geometry.S*Modelica.Constants.g_n));
   Ft[j]*(Re[j]^0.2*geometry.eps^0.6*(1-0.93*0.62)*sin_theta^0.3) = 29.12 * (We[j]*Fr[j])^0.15*geometry.S^0.359;
 end for;

if homotopyMethod.bool_dp and homotopyMethod.useHomotopy then
  // die max-Abfrage braucht man, für den Fall, daß man die Kolonne leer initialisiert und das Inertgas nicht mit abbildet: in dem Fall würde sonst ein Volumenstrom in die falsche Richtung entstehen

  for j in 1:n - 1 loop
 // C[j]=- (p[j] - p[j + 1])/(geometry.H/n)*(1-(0.164 + 71.35*geometry.S)*(4*Ft[j]/geometry.S)^(2/3)  *(3*eta_l_calc[j] * w_sup_l[j]/(rho_l_calc[j]*sin_theta*geometry.eps*9.81*((rho_l_calc[j]-rho_v_calc[j])/rho_l_calc[j]*max(0.1,(1-(p[j] - p[j + 1])/(geometry.H/n)/deltaP_flood)))))^(1/3))^5;
 (p[j] - p[j + 1])=homotopy(actual=- C[j]*(geometry.H/n)/(1-(0.164 + 71.35*geometry.S)*(4*Ft[j]/geometry.S)^(2/3)  *(3*eta_l_calc[j] * w_sup_l[j]/(rho_l_calc[j]*sin_theta*geometry.eps*9.81*((rho_l_calc[j]-rho_v_calc[j])/rho_l_calc[j]*max(0.1,(1-(p[j] - p[j + 1])/(geometry.H/n)/deltaP_flood)))))^(1/3))^5,simplified=homotopyMethod.dp);

 w_sup_v[j] = (-B[j] + sqrt(max(0,B[j]^2 - 4*A[j]*C[j])))/(2*A[j]);
  end for;

//   /*** outlet: half of the n-th discreet element ***/

for j in n:n loop
 // C[j]=- (p[j] - p[j + 1])*2/(geometry.H/n)*(1-(0.164 + 71.35*geometry.S)*(4*Ft[j]/geometry.S)^(2/3)  *(3*eta_l_calc[j] * w_sup_l[j]/(rho_l_calc[j]*sin_theta*geometry.eps*9.81*((rho_l_calc[j]-rho_v_calc[j])/rho_l_calc[j]*max(0.1,(1-(p[j] - p[j + 1])/(geometry.H/n)/deltaP_flood)))))^(1/3))^5;
 (p[j] - p[j + 1])=homotopy(actual=- C[j]/2*(geometry.H/n)/(1-(0.164 + 71.35*geometry.S)*(4*Ft[j]/geometry.S)^(2/3)  *(3*eta_l_calc[j] * w_sup_l[j]/(rho_l_calc[j]*sin_theta*geometry.eps*9.81*((rho_l_calc[j]-rho_v_calc[j])/rho_l_calc[j]*max(0.1,(1-(p[j] - p[j + 1])/(geometry.H/n)/deltaP_flood)))))^(1/3))^5,simplified=homotopyMethod.dp/2);

 w_sup_v[j] = (-B[j] + sqrt(max(0,B[j]^2 - 4*A[j]*C[j])))/(2*A[j]);
end for;

    // A[1]*w_sup_v_in^2 + B[1]*w_sup_v_in - (p_v_in - p[1])*2/(geometry.H/n)*(1-(0.164 + 71.35*geometry.S)*(4*Ft[1]/geometry.S)^(2/3)  *(3*eta_l_calc[1] * w_sup_l[1]/(rho_l_calc[1]*sin_theta*geometry.eps*9.81*((rho_l_calc[1]-rho_v_calc[1])/rho_l_calc[1]*max(0.1,(1-(p_v_in - p[1])*2/(geometry.H/n)/deltaP_flood/2)))))^(1/3))^5=0;
     (p_v_in - p[1])=homotopy(actual=(A[1]*w_sup_v_in^2 + B[1]*w_sup_v_in)/2*(geometry.H/n)/(1-(0.164 + 71.35*geometry.S)*(4*Ft[1]/geometry.S)^(2/3)  *(3*eta_l_calc[1] * w_sup_l[1]/(rho_l_calc[1]*sin_theta*geometry.eps*9.81*((rho_l_calc[1]-rho_v_calc[1])/rho_l_calc[1]*max(0.1,(1-(p_v_in - p[1])*2/(geometry.H/n)/deltaP_flood/2)))))^(1/3))^5,simplified=homotopyMethod.dp/2);

else

  // die max-Abfrage braucht man, für den Fall, daß man die Kolonne leer initialisiert und das Inertgas nicht mit abbildet: in dem Fall würde sonst ein Volumenstrom in die falsche Richtung entstehen

  for j in 1:n - 1 loop
 // C[j]=- (p[j] - p[j + 1])/(geometry.H/n)*(1-(0.164 + 71.35*geometry.S)*(4*Ft[j]/geometry.S)^(2/3)  *(3*eta_l_calc[j] * w_sup_l[j]/(rho_l_calc[j]*sin_theta*geometry.eps*9.81*((rho_l_calc[j]-rho_v_calc[j])/rho_l_calc[j]*max(0.1,(1-(p[j] - p[j + 1])/(geometry.H/n)/deltaP_flood)))))^(1/3))^5;
 (p[j] - p[j + 1])=- C[j]*(geometry.H/n)/(1-(0.164 + 71.35*geometry.S)*(4*Ft[j]/geometry.S)^(2/3)  *(3*eta_l_calc[j] * w_sup_l[j]/(rho_l_calc[j]*sin_theta*geometry.eps*9.81*((rho_l_calc[j]-rho_v_calc[j])/rho_l_calc[j]*max(0.1,(1-(p[j] - p[j + 1])/(geometry.H/n)/deltaP_flood)))))^(1/3))^5;

 w_sup_v[j] = (-B[j] + sqrt(max(0,B[j]^2 - 4*A[j]*C[j])))/(2*A[j]);
  end for;

//   /*** outlet: half of the n-th discreet element ***/

for j in n:n loop
 // C[j]=- (p[j] - p[j + 1])*2/(geometry.H/n)*(1-(0.164 + 71.35*geometry.S)*(4*Ft[j]/geometry.S)^(2/3)  *(3*eta_l_calc[j] * w_sup_l[j]/(rho_l_calc[j]*sin_theta*geometry.eps*9.81*((rho_l_calc[j]-rho_v_calc[j])/rho_l_calc[j]*max(0.1,(1-(p[j] - p[j + 1])/(geometry.H/n)/deltaP_flood)))))^(1/3))^5;
 (p[j] - p[j + 1])=- C[j]/2*(geometry.H/n)/(1-(0.164 + 71.35*geometry.S)*(4*Ft[j]/geometry.S)^(2/3)  *(3*eta_l_calc[j] * w_sup_l[j]/(rho_l_calc[j]*sin_theta*geometry.eps*9.81*((rho_l_calc[j]-rho_v_calc[j])/rho_l_calc[j]*max(0.1,(1-(p[j] - p[j + 1])/(geometry.H/n)/deltaP_flood)))))^(1/3))^5;

 w_sup_v[j] = (-B[j] + sqrt(max(0,B[j]^2 - 4*A[j]*C[j])))/(2*A[j]);
end for;

    // A[1]*w_sup_v_in^2 + B[1]*w_sup_v_in - (p_v_in - p[1])*2/(geometry.H/n)*(1-(0.164 + 71.35*geometry.S)*(4*Ft[1]/geometry.S)^(2/3)  *(3*eta_l_calc[1] * w_sup_l[1]/(rho_l_calc[1]*sin_theta*geometry.eps*9.81*((rho_l_calc[1]-rho_v_calc[1])/rho_l_calc[1]*max(0.1,(1-(p_v_in - p[1])*2/(geometry.H/n)/deltaP_flood/2)))))^(1/3))^5=0;
     (p_v_in - p[1])=(A[1]*w_sup_v_in^2 + B[1]*w_sup_v_in)/2*(geometry.H/n)/(1-(0.164 + 71.35*geometry.S)*(4*Ft[1]/geometry.S)^(2/3)  *(3*eta_l_calc[1] * w_sup_l[1]/(rho_l_calc[1]*sin_theta*geometry.eps*9.81*((rho_l_calc[1]-rho_v_calc[1])/rho_l_calc[1]*max(0.1,(1-(p_v_in - p[1])*2/(geometry.H/n)/deltaP_flood/2)))))^(1/3))^5;
end if;

  annotation (Documentation(info="<html>
</html>", revisions="<html>
<pre><font style=\"color: #006400; \">Documentation last revised: 18.7.2011</font></pre>
</html>"));
end Rocha;
