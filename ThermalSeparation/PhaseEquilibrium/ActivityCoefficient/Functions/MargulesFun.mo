within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions;
function MargulesFun
    input Integer nS=3;
    input Integer k;
    input Modelica.SIunits.MoleFraction x[nS];
    input Real A[nS,nS];
    output Real Margules=0;
protected
   Real margules[nS,nS];

algorithm
    for i in 1:nS loop
     for j in 1:nS loop
      margules [i,j] :=((A[i, k] - 0.5*A[i, j])*x[i]*x[j]);
     end for;
    end for;
    Margules :=sum(margules);
end MargulesFun;
