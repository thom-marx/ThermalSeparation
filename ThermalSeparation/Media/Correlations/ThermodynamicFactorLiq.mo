within ThermalSeparation.Media.Correlations;
package ThermodynamicFactorLiq

  model BaseThermodynamicFactor
    parameter Integer nS= 2;
    constant Real R=Modelica.Constants.R;
    input SI.Temperature T;
    input SI.MoleFraction x[nS];
    output Real Gamma[nS,nS];
  equation

  end BaseThermodynamicFactor;

  model Ideal
    extends BaseThermodynamicFactor;
  equation

      Gamma[:,:]=diagonal(ones(nS));

  end Ideal;

  model NRTL "NRTL"
    extends BaseThermodynamicFactor;
  protected
    Real Q[nS,nS];
    Real epsilon[nS,nS];
    Real tau[nS,nS];
    parameter Real alpha[nS,nS]=fill(2,nS,nS);
    parameter Real g[nS,nS]=fill(2,nS,nS);
    Real G[nS,nS];
    Real S[nS];
    Real C[nS];
  equation

      for i in 1:nS loop
        for m in 1:nS loop
          if i==m then
            tau[i,m]=0;
            G[i,m]=0;
          else
          tau[i,m]=(g[i,m]-g[i,i])/(R*T);
          G[i,m]=exp(-alpha[i,m]*tau[i,m]);
          end if;
        end for;
       end for;

       for i in 1:nS loop
         S[i]=sum(x[k]*G[k,i] for k in 1:nS);
         C[i]=sum(x[k]*G[k,i]*tau[k,i] for k in 1:nS);
       end for;

       for i in 1:nS loop
        for m in 1:nS loop
          epsilon[i,m]=G[i,m]*(tau[i,m]-C[m]/S[m])/S[m];
        end for;
       end for;

       for i in 1:nS loop
        for m in 1:nS loop
          Q[i,m]=epsilon[i,m]+epsilon[m,i]-sum(x[k]*(G[i,k]*epsilon[m,k]+G[m,k]*epsilon[i,k])/S[k] for k in 1:nS);
        end for;
       end for;

       for i in 1:nS loop
        for m in 1:nS loop
          if i==m then
            Gamma[i,m]=1+x[i]*(Q[i,m]-Q[i,nS]);
          else
            Gamma[i,m]=x[i]*(Q[i,m]-Q[i,nS]);
          end if;
         end for;
        end for;

  end NRTL;

  model Test
    parameter Integer nS=5;
    parameter Integer n=3;
    parameter Real x[n,nS]=fill(ones(nS)*1/nS,n);
    //parameter Real R=10;
    parameter Real T[n]={400,300,200};
    //parameter Real V[n,nS]={{10,100},{10,100},{100,10}};
    ThermalSeparation.Media.Correlations.ThermodynamicFactorLiq.NRTL
      nrtlMulticomponent(
      nS=nS,
      n=n,
      x=x,
      T=T);
    ThermalSeparation.Media.Correlations.ThermodynamicFactorLiq.Wilson
      wilsonMulticomponent(                   nS=nS, n=n, x=x, T=T);
    ThermalSeparation.Media.Correlations.ThermodynamicFactorLiq.UNIQUAC
      uniquacMulticomponent(
      nS=nS,
      n=n,
      x=x,
      T=T);
    //UniquacBinary UniquacBinary(nS=nS, n=n,x=x,T=T);
    Ideal ideal(nS=nS, n=n,x=x,T=T);
  equation

  end Test;

  model UNIQUAC "UNIQUAC"
    extends BaseThermodynamicFactor;

    parameter Real q[nS];
    parameter Real r[nS];
    input Real lambda[nS,nS];
  protected
    parameter Integer z=10;
    Real r2;
    Real q2;
    Real Phi[nS];
    Real chi;
    Real Qc[nS,nS];
    Real Qr[nS,nS];
    Real Q[nS,nS];
    Real epsilon[nS,nS];
    Real tau[nS,nS];
    Real S[nS];
  equation

        r2=sum(x[k]*r[k] for k in 1:nS);
        q2=sum(x[k]*q[k] for k in 1:nS);
        chi=sum(x[k] for k in 1:nS);
      for i in 1:nS loop
        Phi[i]=x[i]*q[i]/q2;
        for m in 1:nS loop
          tau[m,i]=exp(-(lambda[m,i]-lambda[m,m])/(R*T));
          Qc[i,m]=-r[i]/r2-r[m]/r2+(r[i]*r[m]/r2^2)*chi-z/2*q2*(r[m]/r2-q[m]/q2)*(r[i]/r2-q[i]/q2);
        end for;
      end for;

       for i in 1:nS loop
         S[i]=sum(Phi[k]*tau[k,i] for k in 1:nS);
         for m in 1:nS loop
           epsilon[i,m]=tau[i,m]/S[m];
           Qr[i,m]=q[i]*q[m]*(1-epsilon[i,m]-epsilon[m,i]+sum(Phi[k]*epsilon[i,k]*epsilon[m,k] for k in 1:nS))/q2;
         end for;
       end for;

       for i in 1:nS loop
         for m in 1:nS loop
           Q[i,m]=Qc[i,m]+Qr[i,m];
         end for;
       end for;

       for i in 1:nS loop
         for m in 1:nS loop
           if i==m then
             Gamma[i,m]=1+x[i]*(Q[i,m]-Q[i,nS]);
           else
             Gamma[i,m]=x[i]*(Q[i,m]-Q[i,nS]);
           end if;
          end for;
        end for;

  end UNIQUAC;

  model Wilson "Wilson"
    extends BaseThermodynamicFactor;

  protected
    parameter Units.DipoleMoment mu[nS,nS]=fill(1,nS,nS);
    parameter Real V[nS]=fill(1,nS);
    Real lambda[nS,nS];
    Real S[nS];
    Real Q[nS,nS];
  equation

    for i in 1:nS loop
      for m in 1:nS loop
        lambda[i,m]=(V[m]/V[i])*exp(-(mu[i,m]-mu[i,i])/(R*T));
      end for;
      S[i]=sum(x[run]*lambda[i,run] for run in 1:nS);
      end for;

    for i in 1:nS loop
      for m in 1:nS loop
        Q[i,m]=-lambda[i,m]/S[i]-lambda[m,i]/S[m]+sum(x[k]*lambda[k,i]*lambda[k,m]/S[k]^2 for k in 1:nS);

      end for;
        end for;

     for i in 1:nS loop
      for m in 1:nS loop
      if i==m then
          Gamma[i,m]=1+x[i]*(Q[i,m]-Q[i,nS]);
      else
          Gamma[i,m]=x[i]*(Q[i,m]-Q[i,nS]);
        end if;
      end for;
      end for;

  end Wilson;
  annotation (Documentation(info="<html>
<p>For the molar flow rate over the phase boundary, the chemical potential &Delta; &mu; is the driving force. In order to relate the driving force to the composition gradient &Delta; x, a correction is needed which is provided by the thermodynamic factors &Gamma;. For the liquid phase, &Gamma; is obtained by the derivation of the activity coefficient with respect to the composition. Taylor and Kooijmann [1] provided this derivation for the Wilson model, the NRTL model and the UNIQUAC model. Those three models are implemented in this package.</p>
<p><br/>References:</p>
<p>[1] R. Taylor, H.A. Kooijmann: Composition derivatives of Activity Models, Chem. Eng. Comm., Vol. 102 (1991), pp. 87-106</p>
</html>"));
end ThermodynamicFactorLiq;
