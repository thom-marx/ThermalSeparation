within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.LiquidPhase;
partial model BaseFugacityCoefficient "base model for fugacity coefficients"
  parameter Integer n=1 annotation(Dialog(enable=false));
  parameter Integer nS=2 
                       annotation(Dialog(enable=false));
  parameter Integer nSL = 2;
  parameter Boolean inertLiquid[nSL]=fill(false,2);
  parameter Integer reorgLiq[nS]={1,2};
replaceable package Medium = 
      ThermalSeparation.Media.WaterBasedLiquid.N2_H2O_O2 
  constrainedby ThermalSeparation.Media.BaseMediumLiquid;

  input SI.Pressure p[n];
  input SI.Temperature T[n];

  output Real phi_sat[n,nSL];//(start=fill(1,n,nS));

protected
  Real phi_sat_aux[n,nS];
equation
  for i in 1:nS loop
    phi_sat[:,reorgLiq[i]] = phi_sat_aux[:,i];
    end for;

  for i in 1:nSL loop
    if inertLiquid[i] then
      phi_sat[:,i] = zeros(n);
    else
    end if;
    end for;

end BaseFugacityCoefficient;
