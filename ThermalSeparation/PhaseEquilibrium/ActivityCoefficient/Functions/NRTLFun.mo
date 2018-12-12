within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions;
function NRTLFun
 input Integer nS;
 input Integer k;
 input Modelica.SIunits.MoleFraction x[nS];
 input Modelica.SIunits.Temperature T;
 input Real alpha[nS,nS];
 input Real g[nS,nS];
 output Real NRTL=0;
protected
 Real nrtl[nS];
 Real tau[nS,nS];
 Real G[nS,nS];

algorithm
   for i in 1:nS loop
    for j in 1:nS loop
     tau[i,j]:=(g[i, j] - g[j, j])/Modelica.Constants.R*T;
     G[i,j]:=if i == j then 1 else -alpha[i, j]*tau[i, j];
    end for;
   end for;
   for i in 1:nS loop
    nrtl[i]:=((x[ i] * G[k, i]/sum(x[ :] .* G[:, i])) * (tau[k, i] - sum(x[ :] .* tau[:, i] .* G[:, i])/sum(x[ :] .* G[:, i])));
   end for;
    NRTL :=        (sum(tau[:, k] .* G[:, k] .* x[ :])/sum(G[:, k] .* x[ :])) + sum(nrtl);
end NRTLFun;
