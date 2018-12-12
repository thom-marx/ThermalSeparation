within ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Liquid;
model ChiltonColburn "Chilton-Colburn analogy between heat and mass transfer"
  extends Liquid.BaseLiquid;

protected
  SI.ThermalConductivity lambda[n] = props.lambda;
  SI.SpecificHeatCapacity cp[n] = props.cp;
equation

  for j in 1:n loop
  //  lambda[j] =Medium.thermalConductivity(state[j]);
  //  cp[j]=Medium.specificHeatCapacityCp(state[j]);

    alpha[j] =max(1000, k_av[j] * rho[j] *cp[j]/(cp[j]*rho[j]/lambda[j]*D_av[j])^(2/3));//if time < 0.1 then 20500 else k_av[j] * rho[j] *cp[j]/(cp[j]*rho[j]/lambda[j]*D_av[j])^(2/3);

    end for;
  annotation (Documentation(info="<html>
<p>[1] Kooijman: Nonequilibrium Column Simulation (Ph.D. Thesis)</p>
</html>"));
end ChiltonColburn;
