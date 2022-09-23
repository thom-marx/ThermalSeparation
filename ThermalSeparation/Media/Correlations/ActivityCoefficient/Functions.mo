within ThermalSeparation.Media.Correlations.ActivityCoefficient;
package Functions
  "Activity coefficient described using functions - often results in index problems!!"

  model BaseActivityCoefficient "base model for activity coefficient"
    parameter Integer n = 1 annotation(Dialog(enable=false));
    parameter Integer nS= 2 annotation(Dialog(enable=false));
    input SI.Temperature T[n];
    input SI.MoleFraction x_l[n,nS];
    output Real gamma[n,nS];
  equation

  end BaseActivityCoefficient;

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
          ThermalSeparation.Media.Correlations.ActivityCoefficient.Functions.MargulesFun(
            nS,
            k,
            x_l[K, :],
            A)/(Modelica.Constants.R*T[K]));
    end for;
    end for;

  end Margules;

  function MargulesFun
      input Integer nS=3;
      input Integer k;
      input Modelica.Units.SI.MoleFraction x[nS];
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

  model NRTL "multicomponent mixture after NRTL eq."
    extends BaseActivityCoefficient;

    parameter Real alpha[nS,nS]={{0.01,0.03},{0.02,0.04}};
    parameter Real g[nS,nS]={{0.01,0.03},{0.04,0.02}};

  equation
   for K in 1:n loop
    for k in 1:nS loop
      gamma[K,k] = exp(NRTLFun(nS,k,x_l[K,:],T[K],alpha,g));
    end for;
   end for;

  end NRTL;

  function NRTLFun
   input Integer nS;
   input Integer k;
   input Modelica.Units.SI.MoleFraction x[nS];
   input Modelica.Units.SI.Temperature T;
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

  model UNIQUAC "multicomponent mixture after UNIQUAC eq."
    extends BaseActivityCoefficient;

    parameter Real q[nS]=fill(1,nS);
    parameter Real r[nS]=fill(3,nS);
    parameter Real u[nS,nS]=fill(2,nS,nS);

  equation
   for K in 1:n loop
    for k in 1:nS loop
        gamma[K, k] =
          ThermalSeparation.Media.Correlations.ActivityCoefficient.Functions.UNIQUACFun(
            nS,
            k,
            x_l[K, :],
            T[K],
            r,
            q,
            u);
    end for;
  end for;

  end UNIQUAC;

  function UNIQUACFun
    parameter Integer z; // z wird hier nicht gesetzt.... sollte 10 sein!
    input Integer nS;
    input Integer k;
    input Modelica.Units.SI.MoleFraction x[nS];
    input Modelica.Units.SI.Temperature T;
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

  model Wilson "multicomponent mixture after Wilson eq."
    extends BaseActivityCoefficient;

    parameter Real Lambda[nS,nS]=fill(1,nS,nS)
      "matrix with the binary coefficients";

  equation
   for K in 1:n loop
    for k in 1:nS loop
        gamma[K, k] =
          ThermalSeparation.Media.Correlations.ActivityCoefficient.Functions.WilsonFun(
            nS,
            k,
            x_l[K, :],
            Lambda);
    end for;
   end for;

  end Wilson;

  function WilsonFun
   input Integer nS;
   input Integer k;
   input Modelica.Units.SI.MoleFraction x[nS];
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
end Functions;
