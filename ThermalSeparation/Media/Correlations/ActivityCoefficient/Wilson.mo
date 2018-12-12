within ThermalSeparation.Media.Correlations.ActivityCoefficient;
model Wilson
  extends BaseActivityCoefficient;

  Real Lambda[nS,nS];
   // "matrix with the binary coefficients";
   parameter Real A[nS,nS]={{0, -696.5, -31.19, 645.7}, {1123.1, 0, 2535.2, 237.5}, {813.2, -547.5, 0, 107.4}, {1917.2, 658, 469.55, 0}};
   parameter Real V_i[nS] = {79.84, 57.54, 44.44, 18.07};

Real term1[nS];
Real term2[nS];
 Real gamma2[nS];
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
  for j in 1:nS loop
    Lambda[i,j] = V_i[j]/V_i[i]*exp(-A[i,j]/(Modelica.Constants.R*T));
  end for;
  end for;
    for i in 1:nS loop
      term1[i] =ln(sum(x_l[k1]*Lambda[i,k1] for k1 in 1:nS));
      term2[i] = sum((x_l[k]*Lambda[k,i])/(sum(x_l[a]*Lambda[k,a] for a in 1:nS)) for k in 1:nS);
      gamma[i]= exp(1-sum((x_l[k]*Lambda[k,i])/(sum(x_l[a]*Lambda[k,a] for a in 1:nS)) for k in 1:nS))/sum(x_l[k1]*Lambda[i,k1] for k1 in 1:nS);
      gamma2[i]=exp(1-term1[i] - term2[i]);
    end for;

end Wilson;
