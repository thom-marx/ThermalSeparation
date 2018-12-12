within ThermalSeparation.HeatAndMassTransfer.ThermodynamicFactorLiq;
model UNIQUACMulticomponent "UNIQUAC"
  extends BaseThermodynamicFactor;
/*
  parameter Integer n = 1;
  parameter Integer nS= 2;
  parameter Real R=10;
  input Real V[n,nS];
  input SI.Temperature T[n];
  input SI.MoleFraction x[n,nS];
  output Real gamma[n,nS-1,nS-1];
*/
protected
  parameter Units.DipoleMoment mu[nS,nS]=fill(1,nS,nS);
  parameter Real q[nS]=fill(1,nS);
  parameter Real r[nS]=fill(1,nS);
  parameter Real lambda[nS,nS]=fill(1,nS,nS);
  parameter Integer z=10;
  Real r2[n];
  Real q2[n];
  Real Phi[n,nS];
  Real chi[n];
  Real Qc[n,nS,nS];
  Real Qr[n,nS,nS];
  Real Q[n,nS,nS];
  Real epsilon[n,nS,nS];
  Real tau[n,nS,nS];
  Real S[n,nS];
equation

  for j in 1:n loop
      r2[j]=sum(x[j,k]*r[k] for k in 1:nS);
      q2[j]=sum(x[j,k]*q[k] for k in 1:nS);
      chi[j]=sum(x[j,k] for k in 1:nS);
    for i in 1:nS loop
      Phi[j,i]=x[j,i]*q[i]/q2[j];
      for m in 1:nS loop
        tau[j,m,i]=exp(-(lambda[m,i]-lambda[m,m])/(R*T[j]));
        Qc[j,i,m]=-r[i]/r2[j]-r[m]/r2[j]+(r[i]*r[m]/r2[j]^2)*chi[j]-z/2*q2[j]*(r[m]/r2[j]-q[m]/q2[j])*(r[i]/r2[j]-q[i]/q2[j]);
      end for;
    end for;
  end for;

   for j in 1:n loop
     for i in 1:nS loop
       S[j,i]=sum(Phi[j,k]*tau[j,k,i] for k in 1:nS);
       for m in 1:nS loop
         epsilon[j,i,m]=tau[j,i,m]/S[j,m];
         Qr[j,i,m]=q[i]*q[m]*(1-epsilon[j,i,m]-epsilon[j,m,i]+sum(Phi[j,k]*epsilon[j,i,k]*epsilon[j,m,k] for k in 1:nS))/q2[j];
       end for;
     end for;
   end for;

   for j in 1:n loop
     for i in 1:nS loop
       for m in 1:nS loop
         Q[j,i,m]=Qc[j,i,m]+Qr[j,i,m];
       end for;
     end for;
   end for;

   for j in 1:n loop
     for i in 1:nS loop
       for m in 1:nS loop
         if i==m then
           Gamma[j,i,m]=1+x[j,i]*(Q[j,i,m]-Q[j,i,nS]);
         else
           Gamma[j,i,m]=x[j,i]*(Q[j,i,m]-Q[j,i,nS]);
         end if;
        end for;
      end for;
   end for;

end UNIQUACMulticomponent;
