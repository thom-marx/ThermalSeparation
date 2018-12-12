within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase;
model PengRobinson
                      //Lüdecke, Lüdecke, Thermodynamik S.566 / Mischungsgrößen: Gmehling & Kolbe Thermodynamik S43.
  extends BaseFugacityCoefficient;
 parameter Real kij[aux2]= fill(0.0267,aux2);
 final parameter Integer aux[:] = {1,3,6,10,15, 21, 28, 36, 45};
 final parameter Integer aux2 = aux[nS-1]
    "number of binary interaction coefficient for Redlich-Kwong equation depending on the number of substances";

    parameter Integer n=1;
    parameter Integer nS=2 "number of components";
    parameter Integer reorgVap[nS] = {1,2};

final parameter SI.Temperature Tcrit[nS]={MediumVapour.Tcrit[reorgVap[i]] for i in 1:nS};
final parameter SI.Pressure pcrit[nS]= {MediumVapour.pcrit[reorgVap[i]] for i in 1:nS};
final parameter Real omega[nS]= {MediumVapour.omega[reorgVap[i]] for i in 1:nS};

protected
  Real b[n]; // Mischungsregel
  Real bi[nS]; // Reinstoffparameter

  Real a[n]; //Gesamt a

  Real aik[n,nS,nS];// Kreuzkoeffizient

  Real alpha[n,nS];

  Real C[n,nS];
  Integer z;
  Real A[nS,nS] "matrix with the binary interaction coefficients";

algorithm
/*** construction of the matrix A using the vector a which is supplied by the user ***/
  for i in 1:nS loop
    for k in i:nS loop
      z:=z+1;
      if i==k then
         A[i,k]:=0;
         z:=z-1;
        else
         A[i,k]:=kij[z];
         A[k,i]:=A[i,k];
      end if;
    end for;
      end for;

equation
// j = Stufe; i/k = Substanz
  for j in 1:n loop
      for i in 1:nS loop
        for k in 1:nS loop
         if i == k then
          aik[j,i,k] = 0.45724*(Modelica.Constants.R^2)*(Tcrit[i]^2)/pcrit[i] * alpha[j,i];
         else

         aik[j,i,k] =(aik[j,i,i]*aik[j,k,k])^(0.5)*(1-A[i,k]);
         end if;
        end for;
      end for;
  end for;

  for j in 1:n loop
    for i in 1:nS loop
     C[j,i] = sum( y[j,k]*aik[j,k,i] for k in 1:nS);
    end for;
   end for;

  for j in 1:n loop
    for i in 1:nS loop
      alpha[j,i] = (1+(0.37464+1.54226*omega[i]-0.26992*omega[i]*omega[i])*(1-(T[j]/Tcrit[i])^(0.5)))^2;
    end for;
  end for;

/*** Mischungsregeln (Van der Waals) ***/

for j in 1:n loop
  a[j] = sum(sum( (y[j,i]*y[j,k]*aik[j,i,k]) for i in 1:nS) for k in 1:nS);
  b[j] = sum((y[j,i]*bi[i]) for i in 1:nS);
end for;

for i in 1:nS loop
  bi[i]=0.0778*(Modelica.Constants.R*Tcrit[i])/pcrit[i];
end for;

  for i in 1:nS loop
    for j in 1:n loop
        phi_aux[j,i] = exp( bi[i]/b[j] * ( p[j]*v[j]/Modelica.Constants.R/T[j] -1) - Modelica.Math.log(p[j]*(v[j]-b[j])/Modelica.Constants.R/T[j]) - a[j]/(2*sqrt(2)*b[j]*Modelica.Constants.R*T[j]) * (2*C[j,i]/a[j] - bi[i]/b[j])*Modelica.Math.log( (v[j] + (1+sqrt(2))*b[j])/(v[j]+(1-sqrt(2))*b[j])));
    end for;
  end for;

end PengRobinson;
