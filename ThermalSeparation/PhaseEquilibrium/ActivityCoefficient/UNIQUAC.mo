within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient;
model UNIQUAC
  extends BaseActivityCoefficient;
  parameter Integer z=10; //Aus Gases and Liquids Seite 257
  parameter Real q[nS]=fill(1,nS);
  parameter Real r[nS]=fill(3,nS);
  parameter Real u[nS,nS]=fill(2,nS,nS);
  input Real delta_u[nS,nS];

protected
 Real phi[nS];
 Real theta[nS];
 Real tau[nS,nS];
 Real l[nS];
 Real term1[nS];
 Real term2[nS];
 Real a;
 Real b;
 Real c;
 Real d;

equation
    a=phi[1];
    b=x_l[1];
    c=theta[1];
    d=q[1];
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
    term2[i]=(1/sum(theta[k]*tau[k,i] for k in 1:nS))^q[i];

    //gamma[j,i]=phi[j,i]/x_l[j,i]*(theta[j,i]/phi[j,i])^(z/2*q[i]);

  end for;

end UNIQUAC;
