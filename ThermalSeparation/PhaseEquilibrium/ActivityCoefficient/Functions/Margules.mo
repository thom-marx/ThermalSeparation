within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions;
model Margules "multicomponent mixture after Margules two-suffix eq."
  extends BaseActivityCoefficient;

  parameter Real a[:] = {100}
    "vector with binary coefficients: for example: {a12, a13, a23}";

protected
  Real A[nS,nS] "matrix with the binary coefficients";
  Integer x;

algorithm
/*** construction of the matrix A using the vector a which is supplied by the user ***/
  for i in 1:nS loop
    for j in i:nS loop
      x:=x+1;
      if i==j then
         A[i,j]:=0;
         x:=x-1;
        else
         A[i,j]:=a[x];
         A[j,i]:=A[i,j];
      end if;
    end for;
  end for;

equation
 for K in 1:n loop
  for k in 1:nS loop
      gamma[K, k] = exp(
        ThermalSeparation.PhaseEquilibrium.ActivityCoefficient.Functions.MargulesFun(
        nS,
        k,
        x_l[K, :],
        A)/(Modelica.Constants.R*T[K]));
  end for;
  end for;

end Margules;
