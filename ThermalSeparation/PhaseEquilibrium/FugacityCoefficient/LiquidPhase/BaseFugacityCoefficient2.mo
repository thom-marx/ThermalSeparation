within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.LiquidPhase;
partial model BaseFugacityCoefficient2 "base model for fugacity coefficients"
  parameter Integer n=1 annotation(Dialog(enable=false));
  parameter Integer nS=2 
                       annotation(Dialog(enable=false));

  input SI.Pressure p;
  input SI.Temperature T;

  output Real phi_sat[nS];//(start=fill(1,n,nS));

end BaseFugacityCoefficient2;
