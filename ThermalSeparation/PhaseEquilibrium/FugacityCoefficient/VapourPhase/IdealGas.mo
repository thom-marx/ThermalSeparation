within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase;
model IdealGas "fugacity coefficient for ideal gas = 1"
  import ThermalSeparation;
  extends BaseFugacityCoefficient;

equation
  phi_aux = ones(n,nS);
end IdealGas;
