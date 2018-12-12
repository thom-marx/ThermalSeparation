within ThermalSeparation.HeatAndMassTransfer.TrayColumn.Vapour;
model Stichlmair
  //extends ThermalSeparation.FilmModel.BaseClasses.TrayColumn.BaseMSVapour;
  extends BaseVapMT;
//Stichlmair: Grundlagen der Dimensionierung des Gas/Flüssigkeit-Kontaktapparates Bodenkolonne, Weinheim, NY: Verlag Chemie, 1978, p. 142
equation

  for j in 1:n loop
    k[j,:] = sqrt(4/Modelica.Constants.pi .* D[j,:] * w_sup[j]/max(1e-5,h[j])/max(1e-5,(1-eps_liq_2ph[j])));

  end for;
end Stichlmair;
