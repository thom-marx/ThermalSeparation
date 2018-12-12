within ThermalSeparation.Media.Correlations.ActivityCoefficient;
model UNIQUAC
  extends BaseActivityCoefficient;
  parameter Integer z=10; //from Gases and Liquids, page 257

//additional inputs to those listed in BaseActivityCoefficient are needed
//this is not a problem, as in a specific medium model one model for the activity coefficient is chosen, which is not supposed to be replaced
  parameter Real q[nS]=fill(1,nS);
  parameter Real r[nS]=fill(3,nS);
  input Real delta_u[nS,nS] "binary interaction energy parameter";

protected
 Real phi[nS];
 Real theta[nS];
 Real tau[nS,nS];
 Real l[nS];
 Real term1[nS];
 Real term2[nS];
equation
for i in 1:nS loop
  l[i]= z/2*(r[i] - q[i]) - (r[i] - 1);
end for;

  for i in 1:nS loop
    theta[i]=q[i]*max(1e-5,x_l[i])/sum(q[k]*x_l[k] for k in 1:nS);
    phi[i]=r[i]*max(1e-5,x_l[i])/sum(r[k]*x_l[k] for k in 1:nS);
    for m in 1:nS loop
      tau[i,m]=exp(-1*(delta_u[i,m])/(Modelica.Constants.R*T));

    end for;

       gamma[i]=min(10.05,   phi[i]/max(1e-5,x_l[i]) * (max(1e-5,theta[i]/phi[i]))^(z/2*q[i]) * term1[i] * term2[i]);

    term1[i]=exp(l[i] + q[i] - phi[i]/max(1e-5,x_l[i])*sum(x_l[k]*l[k] for k in 1:nS) - q[i]*sum(theta[k]*tau[i,k]/sum(theta[a]*tau[a,k] for a in 1:nS) for k in 1:nS));
    term2[i]=(1/max(1e-5,sum(theta[k]*tau[k,i] for k in 1:nS)))^q[i];

    //gamma[j,i]=phi[j,i]/x_l[j,i]*(theta[j,i]/phi[j,i])^(z/2*q[i]);

    end for;

end UNIQUAC;
