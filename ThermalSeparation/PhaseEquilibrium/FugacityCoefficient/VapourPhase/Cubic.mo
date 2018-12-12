within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase;
model Cubic "cubic equation of state, after Redlich-Kwong"
/*** z.B. Reid, Prausnitz: The Properties of Gases & Liquids, 4th edition, McGraw-Hill ***/
  extends BaseFugacityCoefficient;
  /*based on the redlich Kwong constants for cubic equations of state */

  parameter Real kij[nS,nS]= fill(0.01,nS,nS)
    "binary interaction coefficient for Redlich-Kwong equation";

  parameter Integer u=1 "nach Redlich Kwong Näherung";
  parameter Integer w=0 "nach Redlich Kwong Näherung";

protected
  Real h[n,nS] "Hilfsvariable ohne phys. Bedeutung";
  Real h1[n,nS] "quotient aus bi/b";
  Real h2[n] "Hilfsvariable ohne phys. Bedeutung";
  Real delta[n,nS];

  Real Z[n](start=fill(10,n)) "Kompressionsfaktor";
  Real B[n] "Mischkoeffizient";
  Real B1[n] "eigentlich B*";
  Real A1[n] "eigentlich A*";
  Real a[n,nS] "Konstanten nach Redlich Kwong für einzelne Komponenten";
  Real b[n,nS] "Konstanten nach Redlich Kwong für einzelne Komponenten";
  Real am[n] "Konstanten a des Gemischs";
  Real bm[n] "Konstanten b des Gemischs";
Real test[n,nS];
Real test1[n,nS];
Real test2[n,nS];

equation
  for m in 1:n loop
    B[m]=bm[m]-(am[m]/(Modelica.Constants.R*T[m]));
    Z[m]=1+((B[m]*p[m]/1e5)/(Modelica.Constants.R*T[m]));
    B1[m]=(bm[m]*p[m]/1e5)/(Modelica.Constants.R*T[m]);
    A1[m]=(am[m]*p[m]/1e5)/(Modelica.Constants.R^2*T[m]^2);
    am[m]=sum(sum((y[m,i]*y[m,j])*(a[m,i]*a[m,j])^(0.5)*(1-kij[i,j]) for i in 1:nS) for j in 1:nS);
    bm[m]=sum(y[m,i]*b[m,i] for i in 1:nS);
    h2[m]=sum(sum(x[m,j]*a[m,j]^(0.5)*(1-kij[k,j]) for j in 1:nS) for k in 1:nS);

    for i in 1:nS loop
      a[m,i]=  (0.42748*Modelica.Constants.R^2*MediumVapour.Tcrit[i]^(2.5))/(T[m]^(0.5)*
          MediumVapour.pcrit[i]/1e5);
      b[m,i]=  0.08664*Modelica.Constants.R*MediumVapour.Tcrit[i]/(MediumVapour.pcrit[
          i]/1e5);

      h[m,i]=sum((y[m, j]*MediumVapour.Tcrit[j])/(MediumVapour.pcrit[j]/1e5)         for j in 1:nS);

      h1[m,i]=(MediumVapour.Tcrit[i]/(MediumVapour.pcrit[i]/1e5))        /h[m,i];

      delta[m,i]=((2*a[m,i]^(0.5))/am[m])*h2[m];
test[m,i] =  exp(h1[m,i]*(Z[m] - 1) - Modelica.Math.log(max(1e-5,(Z[m] - B1[m]))) + (A1[m]/(B1[m]*sqrt(u^2
         - 4*w)))*(h1[m,i] - delta[m,i])*Modelica.Math.log(max(1e-5,((2*Z[m] + B1[m]*(u + sqrt(u^2 - 4*w)))/
        (2*Z[m] + B1[m]*(u - sqrt(u^2 - 4*w)))))));
test1[m,i] = Z[m] - B1[m];
test2[m,i] = (2*Z[m] + B1[m]*(u + sqrt(u^2 - 4*w)))/(2*Z[m] + B1[m]*(u - sqrt(u^2 - 4*w)));
      phi_aux[m,i]=exp(h1[m,i]*(Z[m] - 1) - Modelica.Math.log(max(1e-5,(Z[m] - B1[m]))) + (A1[m]/(B1[m]*sqrt(u^2
         - 4*w)))*(h1[m,i] - delta[m,i])*Modelica.Math.log(max(1e-5,((2*Z[m] + B1[m]*(u + sqrt(u^2 - 4*w)))/
        (2*Z[m] + B1[m]*(u - sqrt(u^2 - 4*w)))))));
    end for;
  end for;

end Cubic;
