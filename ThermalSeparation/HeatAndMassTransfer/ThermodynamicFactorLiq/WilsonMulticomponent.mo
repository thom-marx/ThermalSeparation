within ThermalSeparation.HeatAndMassTransfer.ThermodynamicFactorLiq;
model WilsonMulticomponent "Wilson"
  extends BaseThermodynamicFactor;

protected
  parameter Units.DipoleMoment mu[nS,nS]=fill(1,nS,nS);
  parameter Real V[n,nS]=fill(1,n,nS);
  Real lambda[n,nS,nS];
  Real S[n,nS];
  Real Q[n,nS,nS];
equation
  for j in 1:n loop
  for i in 1:nS loop
    for m in 1:nS loop
      lambda[j,i,m]=(V[j,m]/V[j,i])*exp(-(mu[i,m]-mu[i,i])/(R*T[j]));
    end for;
    S[j,i]=sum(x[j,run]*lambda[j,i,run] for run in 1:nS);
    end for;

  for i in 1:nS loop
    for m in 1:nS loop
      Q[j,i,m]=-lambda[j,i,m]/S[j,i]-lambda[j,m,i]/S[j,m]+sum(x[j,k]*lambda[j,k,i]*lambda[j,k,m]/S[j,k]^2 for k in 1:nS);

    end for;
      end for;

   for i in 1:nS loop
    for m in 1:nS loop
    if i==m then
        Gamma[j,i,m]=1+x[j,i]*(Q[j,i,m]-Q[j,i,nS]);
    else
        Gamma[j,i,m]=x[j,i]*(Q[j,i,m]-Q[j,i,nS]);
      end if;
    end for;
    end for;
    end for;

end WilsonMulticomponent;
