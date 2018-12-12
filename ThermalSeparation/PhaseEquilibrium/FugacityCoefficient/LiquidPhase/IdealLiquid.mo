within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.LiquidPhase;
model IdealLiquid "fugacity coefficient for ideal liquid = 1"
  import ThermalSeparation;
  extends
    ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase.BaseFugacityCoefficient;
equation
  phi_sat_aux = ones(n,nS);
end IdealLiquid;
