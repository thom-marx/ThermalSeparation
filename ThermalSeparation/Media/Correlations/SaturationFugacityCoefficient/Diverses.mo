within ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient;
package Diverses
                 //Unfertige Sachen, vollendete in Vaporphase

  model PengRobinson
    extends ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient.BaseFugacityCoefficientxxx;
  parameter Real kij[nS,nS]= fill(0.01,nS,nS)
      "binary interaction coefficient for Redlich-Kwong equation";

  protected
    Real b; // Mischungsregel
    Real bk[nS]; // Reinstoffparameter

    Real a[n]; //Gesamt a

    Real aij[n,nS,nS];// Kreuzkoeffizient

    Real am[n,nS];  // Reinstoffparameter
    Real alpha[n,nS];

    Real C[n,nS];
  equation
  // m = Stufe; k = Substanz

  for m in 1:n loop
    for k in 1:nS loop
     C[m,k] = sum( x[m,f]*aij[m,f,k] for f in 1:nS);
    end for;
  end for;

    for m in 1:n loop
      for i in 1:nS loop
        for j in 1:nS loop
         aij[m,i,j] =(am[m,i]*am[m,j])^(0.5)*(1-kij[i,j]);
        end for;
      end for;
    end for;

    for k in 1:nS loop
     for m in 1:n loop
      am[m,k] = 0.45724*(Modelica.Constants.R^2)*(MediumLiquid.Tcrit[k]^2)/MediumLiquid.pcrit[k] * alpha[m,k];
     end for;
    end for;

    for k in 1:nS loop
      for m in 1:n loop
        alpha[m,k] = (1+(0.37464+1.54226*MediumLiquid.omega[k]-0.26992*MediumLiquid.omega[k]*MediumLiquid.omega[k])*(1-(T[m]/MediumLiquid.Tcrit[k])^(0.5)))^2;
      end for;
    end for;

  /*** Mischungsregeln ***/
    b = sum(sum( (x[m,k]*bk[k]) for k in 1:nS) for m in 1:n);
  for m in 1:n loop
    a[m] = sum(sum( (x[m,i]*x[m,j]*aij[m,i,j]) for i in 1:nS) for j in 1:nS);
  end for;

  for k in 1:nS loop
    bk[k]=0.0778*(Modelica.Constants.R*MediumLiquid.Tcrit[k])/MediumLiquid.pcrit[k];
  end for;

    for k in 1:nS loop
      for m in 1:n loop

          phi[m,k] = exp(-Modelica.Math.log(p[m]/1e5*v[m]/(Modelica.Constants.R*T[m])-p[m]/1e5*b/(Modelica.Constants.R*T[m])) + (p[m]/1e5*v[m]/(Modelica.Constants.R*T[m])-1)*bk[k]/b - (a[m]/(2^(1.5)*b*Modelica.Constants.R*T[m]))*(1/a[m]*(2*C[m,k])-bk[k]/b)*Modelica.Math.log((v[m]+(2^(0.5)+1)*b)/(v[m]-(2^(0.5)-1)*b)));
      end for;
    end for;

  end PengRobinson;

  model SoaveRedlichKwong
                          // Buch: Thermodynamik (Gmehling, Kolbe), S. 97 equ. 3.75, S. 43
    extends ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient.BaseFugacityCoefficientxxx;

  parameter Real kij[nS,nS]= fill(0.01,nS,nS)
      "binary interaction coefficient for Redlich-Kwong equation";

  protected
    Real b; // Mischungsregel
    Real bk[nS]; // Reinstoffparameter

    Real a[n]; //Gesamt a

    Real aij[n,nS,nS];// Kreuzkoeffizient

    Real am[n,nS];  // Reinstoffparameter
    Real alpha[n,nS];

    Real C[n,nS];
  equation
  // m = Stufe; k = Substanz

  for m in 1:n loop
    for k in 1:nS loop
     C[m,k] = sum( x[m,f]*aij[m,f,k] for f in 1:nS);
    end for;
  end for;

    for m in 1:n loop
      for i in 1:nS loop
        for j in 1:nS loop
         aij[m,i,j] =(am[m,i]*am[m,j])^(0.5)*(1-kij[i,j]);
        end for;
      end for;
    end for;

    for k in 1:nS loop
     for m in 1:n loop
      am[m,k] = 0.42748*(Modelica.Constants.R^2)*(MediumLiquid.Tcrit[k]^2)/MediumLiquid.pcrit[k] * alpha[m,k];
     end for;
    end for;

    for k in 1:nS loop
      for m in 1:n loop
        alpha[m,k] = (1+(0.48+1.574*MediumLiquid.omega[k]-0.176*MediumLiquid.omega[k]*MediumLiquid.omega[k])*(1-(T[m]/MediumLiquid.Tcrit[k])^(0.5)))^2;
      end for;
    end for;

  /*** Mischungsregeln ***/
    b = sum(sum( (x[m,k]*bk[k]) for k in 1:nS) for m in 1:n);
  for m in 1:n loop
    a[m] = sum(sum( (x[m,i]*x[m,j]*aij[m,i,j]) for i in 1:nS) for j in 1:nS);
  end for;

  for k in 1:nS loop
    bk[k]=0.08664*(Modelica.Constants.R*MediumLiquid.Tcrit[k])/MediumLiquid.pcrit[k];
  end for;

    for k in 1:nS loop
      for m in 1:n loop

          phi[m,k] = exp( Modelica.Math.log(v[m]/(v[m]-b)) - (2*C[m,k]/(Modelica.Constants.R*T[n]*b))*Modelica.Math.log((v[n]+b)/v[n]) + bk[k]/(v[n]-b) - Modelica.Math.log((p[n]/1e5)*v[n]/(Modelica.Constants.R*T[n])) + (a[m]*bk[k]/(Modelica.Constants.R*T[n]*b^2))*(Modelica.Math.log((v[n]+b)/v[n] - b/(v[n]+b))));
      end for;
    end for;

  end SoaveRedlichKwong;

  model SoaveRedlichKwong2
                          // Buch: Thermodynamik (Gmehling, Kolbe), S. 97 equ. 3.75, S. 43
    extends ThermalSeparation.Media.Correlations.SaturationFugacityCoefficient.BaseFugacityCoefficientxxx;

  parameter Real kij[nS,nS]= fill(0.01,nS,nS)
      "binary interaction coefficient for Redlich-Kwong equation";

  protected
    Real b; // Mischungsregel
    Real bk[nS]; // Reinstoffparameter

    Real a[n]; //Gesamt a

    Real aij[n,nS,nS];// Kreuzkoeffizient

    Real am[n,nS];  // Reinstoffparameter
    Real alpha[n,nS];

    Real C[n,nS];
  equation
  // m = Stufe; k = Substanz

  for m in 1:n loop
    for k in 1:nS loop
     C[m,k] = sum( x[m,f]*aij[m,f,k] for f in 1:nS);
    end for;
  end for;

    for m in 1:n loop
      for i in 1:nS loop
        for j in 1:nS loop
         aij[m,i,j] =(am[m,i]*am[m,j])^(0.5)*(1-kij[i,j]);
        end for;
      end for;
    end for;

    for k in 1:nS loop
     for m in 1:n loop
      am[m,k] = 0.42748*(Modelica.Constants.R^2)*(MediumLiquid.Tcrit[k]^2)/MediumLiquid.pcrit[k] * alpha[m,k];
     end for;
    end for;

    for k in 1:nS loop
      for m in 1:n loop
        alpha[m,k] = (1+(0.48+1.574*MediumLiquid.omega[k]-0.176*MediumLiquid.omega[k]*MediumLiquid.omega[k])*(1-(T[m]/MediumLiquid.Tcrit[k])^(0.5)))^2;
      end for;
    end for;

  /*** Mischungsregeln ***/
    b = sum(sum( (x[m,k]*bk[k]) for k in 1:nS) for m in 1:n);
  for m in 1:n loop
    a[m] = sum(sum( (x[m,i]*x[m,j]*aij[m,i,j]) for i in 1:nS) for j in 1:nS);
  end for;

  for k in 1:nS loop
    bk[k]=0.08664*(Modelica.Constants.R*MediumLiquid.Tcrit[k])/MediumLiquid.pcrit[k];
  end for;

    for k in 1:nS loop
      for m in 1:n loop

          phi[m,k] = exp(-Modelica.Math.log(p[m]/1e5*v[m]/(Modelica.Constants.R*T[m])-p[m]/1e5*b/(Modelica.Constants.R*T[m])) + (p[m]/1e5*v[m]/(Modelica.Constants.R*T[m])-1)*bk[k]/b - (a[m]/(b*Modelica.Constants.R*T[m]))*(1/a[m]*(2*C[m,k])-bk[k]/b)*Modelica.Math.log(1+b/v[m]));
      end for;
    end for;

  end SoaveRedlichKwong2;
end Diverses;
