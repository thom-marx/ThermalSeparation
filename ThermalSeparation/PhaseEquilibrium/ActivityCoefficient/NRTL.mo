within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient;
model NRTL
  extends BaseActivityCoefficient;

  parameter Real alpha[nS,nS]=fill(0.1,nS,nS);//{{0.01,0.03},{0.02,0.04}};
  parameter Real g[nS,nS]=fill(0.01,nS,nS);//{{0.01,0.03},{0.04,0.02}};

protected
  Real term1[nS];
  Real term2[nS];
  Real tau[nS,nS];
  Real G[nS,nS];
equation

   for i in 1:nS loop
    for m in 1:nS loop
     tau[i,m]=(g[i, m] - g[m, m])/Modelica.Constants.R*T;
     G[i,m]=if i == m then 1 else -alpha[i, m]*tau[i, m];
    end for;
   end for;

  for i in 1:nS loop
      term1[i]=sum(tau[k,i]*G[k,i]*x_l[k] for k in 1:nS)/sum(G[k,i]*x_l[k] for k in 1:nS);
      term2[i]=sum((x_l[k]*G[i,k]/sum(G[k1,k]*x_l[k] for k1 in 1:nS))*(tau[i,k]-sum(x_l[k1]*tau[k1,k]*G[k1,k] for k1 in 1:nS)/sum(x_l[k1]*G[k1,k] for k1 in 1:nS)) for k in 1:nS);
      gamma[i]=exp(term1[i]+term2[i]);
  end for;

end NRTL;
