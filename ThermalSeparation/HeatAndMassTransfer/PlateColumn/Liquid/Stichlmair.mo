within ThermalSeparation.HeatAndMassTransfer.PlateColumn.Liquid;
model Stichlmair
  extends Base;
//Stichlmair: Grundlagen der Dimensionierung des Gas/Flüssigkeit-Kontaktapparates Bodenkolonne, Weinheim, NY: Verlag Chemie, 1978, p. 142

equation
    // Werte belegen, damit das System nicht singulär wird:
//   Re[:]=fill(0,n);
//   Sc[:,:]=fill(0,n,a);
//   Sh[:,:]=fill(0,n,a);
  for j in 1:n loop
    for i in 1:a loop
   // k_2[j,:] = sqrt(4/Modelica.Constants.pi .* D[j,:] * w_sup[j]/max(1e-5,abs(h[j]))/max(1e-5,abs((1-eps_liq_2ph))));
    k_2[j,i] = max(1e-5,sqrt(4/Modelica.Constants.pi .* D[j,i] * w_sup[j]/max(1e-5,abs(h[j]))/max(1e-5,abs((1-eps_liq_2ph[j])))));
    //k_2[j,:] = sqrt(4/Modelica.Constants.pi .* D[j,:] * w_sup[j]);
    end for;
  end for;

                                                  annotation(choicesAllMatching);
end Stichlmair;
