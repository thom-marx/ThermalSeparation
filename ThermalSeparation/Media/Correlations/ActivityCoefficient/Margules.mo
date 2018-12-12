within ThermalSeparation.Media.Correlations.ActivityCoefficient;
model Margules
  extends BaseActivityCoefficient;
    parameter Integer aux[:] = {1,3,6,10,15, 21, 28, 36, 45};
  parameter Integer aux2 = aux[nS-1]
    "number of binary mass transfer coefficients depending on the number of substances";

 // parameter Real a[nS+5] = fill(100,nS+5)
  //  "vector with binary coefficients: for example: {a12, a13, a23}";
  //                                        //hier nochmal gucken.. wusste nicht genau was das ist und hab es so geändert dass es läuft!
  parameter Real a[aux2] = fill(100,aux2);
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

   for i in 1:nS loop
     gamma[i]=exp(sum(sum((A[a,i]-0.5*A[a,b])*x_l[a]*x_l[b] for b in 1:nS) for a in 1:nS)/(Modelica.Constants.R*T));
   end for;

end Margules;
