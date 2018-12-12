within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions;
model NRTL "multicomponent mixture after NRTL eq."
  extends BaseActivityCoefficient;

  parameter Real alpha[nS,nS]={{0.01,0.03},{0.02,0.04}};
  parameter Real g[nS,nS]={{0.01,0.03},{0.04,0.02}};

equation
 for K in 1:n loop
  for k in 1:nS loop
    gamma[K,k] = exp(NRTLFun(nS,k,x_l[K,:],T[K],alpha,g));
  end for;
 end for;

end NRTL;
