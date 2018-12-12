within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase;
model SoaveRedlichKwong
                        // Buch: Thermodynamik (Gmehling, Kolbe), S. 97 equ. 3.75, S. 43
  extends BaseFugacityCoefficient;

 final parameter Integer aux[:] = {1,3,6,10,15, 21, 28, 36, 45};
 final parameter Integer aux2 = aux[nS-1]
    "number of binary interaction coefficient for Redlich-Kwong equation depending on the number of substances";

  parameter Real kij[aux2]= fill(0.0267,aux2)
    "binary interaction coefficient for Redlich-Kwong equation";

    parameter Integer n=1;
    parameter Integer nS=2 "number of components";
    parameter Integer reorgVap[nS] = {1,2};

final parameter SI.Temperature Tcrit[nS]={MediumVapour.Tcrit[reorgVap[i]] for i in 1:nS};
final parameter SI.Pressure pcrit[nS]= {MediumVapour.pcrit[reorgVap[i]] for i in 1:nS};
final parameter Real omega[nS]= {MediumVapour.omega[reorgVap[i]] for i in 1:nS};
protected
  Real b[n]; // Mischungsregel  [m3/mol]
  Real bi[nS]; // Reinstoffparameter

  Real a[n]; //Gesamt a

  Real aik[n,nS,nS];// Kreuzkoeffizient

    // Reinstoffparameter [N*m^4/mol^2]
  Real alpha[n,nS];

  Real C[n,nS];

  Integer z;
 Real A[nS,nS] "matrix with the binary interaction coefficients";

algorithm
/*** construction of the matrix A using the vector kij which is supplied by the user ***/
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

// Druck in Pa und molares Volumen in m3/mol
for j in 1:n loop
    for i in 1:nS loop
      for k in 1:nS loop
       if i == k then
        aik[j,i,k] = 0.42748*((Modelica.Constants.R^2)*(Tcrit[i]^2))/(pcrit[i]) * alpha[j,i];
       else

       aik[j,i,k] =(aik[j,i,i]*aik[j,k,k])^(0.5)*(1-A[i,k]);
       end if;
      end for;
    end for;
end for;

/*** Mischungsregeln ***/

for j in 1:n loop
  a[j] = sum(sum( (y[j,i]*y[j,k]*aik[j,i,k]) for i in 1:nS) for k in 1:nS);
end for;

for i in 1:nS loop
  bi[i]=0.08664*(Modelica.Constants.R*Tcrit[i])/(pcrit[i]);
end for;

for j in 1:n loop
b[j] =sum((y[j,i]*bi[i]) for i in 1:nS);
end for;

  for j in 1:n loop
    for i in 1:nS loop

       alpha[j,i] = (1+(0.48+1.574*omega[i]-0.176*omega[i]^2)*(1-((T[j]/Tcrit[i])^(0.5))))^2;

       C[j,i] = sum( y[j,k]*aik[j,k,i] for k in 1:nS);

    phi_aux[j,i] =  exp( Modelica.Math.log(v[j]/(v[j]-b[j])) - (2*C[j,i]/(Modelica.Constants.R*T[j]*b[j]))*Modelica.Math.log((v[j]+b[j])/v[j]) + bi[i]/(v[j]-b[j]) - Modelica.Math.log((p[j])*v[j]/(Modelica.Constants.R*T[j])) + (a[j]*bi[i]/(Modelica.Constants.R*T[j]*b[j]^2))*(Modelica.Math.log((v[j]+b[j])/v[j]) - b[j]/(v[j]+b[j])));

    end for;
  end for;

end SoaveRedlichKwong;
