within ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Vapour;
model ChiltonColburn "Chilton-Colburn analogy between heat and mass transfer"
  extends BaseVapour;
 // ThermalSeparation.Units.CoefficentOfMassTransfer k_av[n]
   // "averaged mass transfer coefficient";
 // SI.DiffusionCoefficient D_av[n] "averaged diffusion coefficient";
protected
  SI.ThermalConductivity lambda[n] = props.lambda;
  SI.SpecificHeatCapacity cp[n] = props.cp;
equation

  for j in 1:n loop
   // lambda[j] = Medium.thermalConductivity(state[j]);
  //  cp[j]=Medium.specificHeatCapacityCp(state[j]);
    alpha[j] =  max(1,k_av[j] * rho[j] *cp[j]/(cp[j]*rho[j]/lambda[j]*D_av[j])^(2/3));//if time < 0.1 then 300 else k_av[j] * rho[j] *cp[j]/(cp[j]*rho[j]/lambda[j]*D_av[j])^(2/3);
    end for;
  annotation (Documentation(info="<html>
<p>[1] Kooijman: Nonequilibrium Column Simulation (Ph.D. Thesis)</p>
</html>"));
end ChiltonColburn;
