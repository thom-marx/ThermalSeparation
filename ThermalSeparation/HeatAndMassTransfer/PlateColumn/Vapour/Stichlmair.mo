within ThermalSeparation.HeatAndMassTransfer.PlateColumn.Vapour;
model Stichlmair
  extends Base;
//Stichlmair: Grundlagen der Dimensionierung des Gas/Flüssigkeit-Kontaktapparates Bodenkolonne, Weinheim, NY: Verlag Chemie, 1978, p. 142
equation
    // Werte belegen, damit das System nicht singulär wird:
//   Re[:]=fill(0,n);
//   Sc[:,:]=fill(0,n,a);
//   Sh[:,:]=fill(0,n,a);
  for j in 1:n loop
    k_2[j,:] = sqrt(4/Modelica.Constants.pi .* D[j,:] * w_sup[j]/max(1e-5,h[j])/max(1e-5,(1-eps_liq_2ph[j])));
  end for;
end Stichlmair;
