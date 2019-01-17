within ThermalSeparation.PressureLoss.SprayColumn;
model PipeFlow "calculation of Vdot using the particle model and zeta values"

  extends BasicPressureLossSpray;

  SI.Density rho_v_av[n-1];

  parameter SI.Height h = geometry.H/n "height of one discreet element";
  parameter Real lambda_const = 10
    "lambda value must be chosen higher than in reality, otherwise the pressure loss is far too small and the simulation fails";

protected
  Real lambda[n-1];
  Real lambda_in;
  Real lambda_out;
  // Real Re[n-1];
  // Real Re_in;
  // Real Re_out;

equation
//Gleichung aus Engel: Fluiddynamik in Füllkörper- und Packungskolonnen ..., Gl. (10); in der Gleichung (3.13b) von Forner ist wohl ein Fehler (h_j)

// die max-Abfrage braucht man, für den Fall, daß man die Kolonne leer initialisiert und das Inertgas nicht mit abbildet: in dem Fall würde sonst ein Volumenstrom in die falsche Richtung entstehen

for j in 1:n-1 loop
  rho_v_av[j] = (rho_v_calc[j]+rho_v_calc[j+1])/2;

lambda[j] =lambda_const;// 0.3164/Re[j];
end for;

for j in 1:n-1 loop
Vdot[j] = sign(p[j]-p[j+1])*sqrt(abs(p[j]-p[j+1]) *2 /rho_v_av[j] /h*geometry.r1*2/lambda[j]*geometry.A^2);
end for;

/*** outlet: half of the n-th discreet element ***/

lambda_out =lambda_const;//  0.3164/Re_out;
Vdot[n] =  sign(p[n]-p[n+1])*sqrt(max(1e-5,abs(p[n]-p[n+1]))*2 /rho_v_calc[n] /(h/2)*geometry.r1*2/lambda_out*geometry.A^2);

/*** entry: half of the first discreet element ***/
lambda_in =lambda_const;//  0.3164/Re_in;
p_v_in = p[1]+ sign(Vdot_in)*Vdot_in*Vdot_in/geometry.A^2 * rho_v_calc[1]/2 * lambda_in * (h/2)/geometry.r1/2;//max(1e-10,f_in/ (6*hu_dyn[1] / d_L_in + geometry.a)/ rho_v_calc[1] * (geometry.A_free)^2))*(h/2);

assert(not homotopyMethod.bool_dp,"this pressure loss model does not support homotopy - please deactivate");

end PipeFlow;
