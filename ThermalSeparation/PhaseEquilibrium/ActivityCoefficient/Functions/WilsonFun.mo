within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions;
function WilsonFun
 input Integer nS;
 input Integer k;
 input Modelica.SIunits.MoleFraction x[nS];
 input Real Lambda[nS,nS];

 output Real Wilson=0;

protected
 Real wilson[nS];

algorithm
 for i in 1:nS loop
  wilson[i]:=x[i]*Lambda[i, k]/sum(x[:] .* Lambda[i, :]);
 end for;
 Wilson:=exp(1 - sum(wilson))/(sum(x[:] .* Lambda[k, :]));

end WilsonFun;
