within ThermalSeparation.Media.Correlations.ActivityCoefficient;
model NRTL
  extends BaseActivityCoefficient;

//additional inputs to those listed in BaseActivityCoefficient are needed
//this is not a problem, as in a specific medium model one model for the activity coefficient is chosen, which is not supposed to be replaced
 input Real alpha[nS,nS] "non-randomness parameters";
 input Real g[nS,nS]=fill(0.01,nS,nS);

protected
  Real term1[nS];
  Real term2[nS];
  Real tau[nS,nS];
  Real G[nS,nS];
  Real test[nS];
equation

   for i in 1:nS loop
    for m in 1:nS loop
     tau[i,m]=g[i, m]/(Modelica.Constants.R*T);
     // tau[i,m]=(g[i, m] - g[m, m])/Modelica.Constants.R*T;
     G[i,m]=if i == m then 1 else exp(-alpha[i, m]*tau[i, m]);
    end for;
   end for;

  for i in 1:nS loop
    test[i]=max(1e-6,sum(G[k1,i]*x_l[i] for k1 in 1:nS));
      term1[i]=sum(tau[k,i]*G[k,i]*x_l[k] for k in 1:nS)/sum(G[k,i]*x_l[k] for k in 1:nS);
      term2[i]=sum((x_l[k]*G[i,k]/test[k]) * (tau[i,k]-sum(x_l[k1]*tau[k1,k]*G[k1,k] for k1 in 1:nS)/sum(x_l[k1]*G[k1,k] for k1 in 1:nS)) for k in 1:nS);
      // term2[i]=sum((x_l[k]*G[i,k]/sum(G[k1,k]*x_l[k] for k1 in 1:nS))*(tau[i,k]-sum(x_l[k1]*tau[k1,k]*G[k1,k] for k1 in 1:nS)/sum(x_l[k1]*G[k1,k] for k1 in 1:nS)) for k in 1:nS);
      gamma[i]=exp(term1[i]+term2[i]);
  end for;

end NRTL;
