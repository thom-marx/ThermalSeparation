within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions;
function UNIQUACFun
  parameter Integer z; // z wird hier nicht gesetzt.... sollte 10 sein!
  input Integer nS;
  input Integer k;
  input Modelica.SIunits.MoleFraction x[nS];
  input Modelica.SIunits.Temperature T;
  input Real r[nS];
  input Real q[nS];
  input Real u[nS,nS];
  output Real UNIQUAC;

protected
 Real phi[nS];
 Real theta[nS];
 Real tau[nS,nS];
 Real l[nS];
 Real uniquac[nS];

algorithm
 for i in 1:nS loop
  for j in 1:nS loop
   tau[i,j]:=exp(-1*(u[i, j] - u[j, j])/Modelica.Constants.R*T);// OK
  end for;
  l[i]:=z/2*(r[i] - q[i]) - (r[i] - 1);//OK
  phi[i]:=q[i]*x[i]/sum(q[:] .* x[:]);//OK
  theta[i]:=r[i]*x[i]/sum(r[:] .* x[:]);//OK
  uniquac[i]:=phi[i]*tau[i, k]/sum(phi[:] .* tau[:, i]);
 end for;
 UNIQUAC :=theta[k]/x[k]*(phi[k]/theta[k])^(z/2*q[k])*exp(l[k] + q[k] - theta[k]/x[k]*sum(x[:] .* l[:]) - q[k]*sum(uniquac))/sum(phi[:].*tau[:,k])^q[k];
end UNIQUACFun;
