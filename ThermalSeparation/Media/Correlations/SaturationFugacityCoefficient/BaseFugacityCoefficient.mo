within ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient;
partial model BaseFugacityCoefficient "base model for fugacity coefficients"
  parameter Integer n=1 annotation(Dialog(enable=false));
  parameter Integer nS=2
                       annotation(Dialog(enable=false));

  parameter Real omega[nS];
  parameter SI.Temperature Tcrit[nS];
  parameter SI.Pressure pcrit[nS];
  parameter Real mu[nS];
  parameter Integer NoOfEq[nS];

  input SI.Pressure p_sat[nS];
  input SI.Pressure p;
  input SI.Temperature T;

  output Real phi_sat[nS];//(start=fill(1,n,nS));

end BaseFugacityCoefficient;
