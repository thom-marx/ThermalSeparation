within ThermalSeparation.HeatAndMassTransfer.ThermodynamicFactorLiq;
model NRTLMulticomponent "NRTL"
  extends BaseThermodynamicFactor;
protected
  Real Q[n,nS,nS];
  Real epsilon[n,nS,nS];
  Real tau[n,nS,nS];
  parameter Real alpha[nS,nS]=fill(2,nS,nS);//{{2,2},{2,2}};
  parameter Real g[nS,nS]=fill(2,nS,nS);//{{2,2},{2,2}};
  Real G[n,nS,nS];
  Real S[n,nS];
  Real C[n,nS];
equation
  for j in 1:n loop
    for i in 1:nS loop
      for m in 1:nS loop
        if i==m then
          tau[j,i,m]=0;
          G[j,i,m]=0;
        else
        tau[j,i,m]=(g[i,m]-g[i,i])/(R*T[n]);
        G[j,i,m]=exp(-alpha[i,m]*tau[j,i,m]);
        end if;
      end for;
     end for;

     for i in 1:nS loop
       S[j,i]=sum(x[j,k]*G[j,k,i] for k in 1:nS);
       C[j,i]=sum(x[j,k]*G[j,k,i]*tau[j,k,i] for k in 1:nS);
     end for;

     for i in 1:nS loop
      for m in 1:nS loop
        epsilon[j,i,m]=G[j,i,m]*(tau[j,i,m]-C[j,m]/S[j,m])/S[j,m];
      end for;
     end for;

     for i in 1:nS loop
      for m in 1:nS loop
        Q[j,i,m]=epsilon[j,i,m]+epsilon[j,m,i]-sum(x[j,k]*(G[j,i,k]*epsilon[j,m,k]+G[j,m,k]*epsilon[j,i,k])/S[j,k] for k in 1:nS);
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

end NRTLMulticomponent;
