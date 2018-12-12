within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient;
model Wilson
  extends BaseActivityCoefficient;

  parameter Real Lambda[nS,nS]=fill(1.1,nS,nS)
    "matrix with the binary coefficients";
 //Real wilson[n,nS,nS];
// Real gamma2[n,nS];
equation
/*  
 for K in 1:n loop
  for k in 1:nS loop
    for i in 1:nS loop
       wilson[K,k,i]=x_l[K,i]*Lambda[i, k]/sum(x_l[K,:] .* Lambda[i, :]);
    end for;
    gamma[K,k] =min(1.05, exp(1 - sum(wilson[K,k,:]))/(sum(x_l[K,:] .* Lambda[k, :])));
  end for;
 end for;
 */

    for i in 1:nS loop
      gamma[i]= exp(1-sum((x_l[k]*Lambda[k,i])/(sum(x_l[a]*Lambda[k,a] for a in 1:nS)) for k in 1:nS))/sum(x_l[k1]*Lambda[i,k1] for k1 in 1:nS);
    end for;

end Wilson;
