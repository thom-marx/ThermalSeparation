within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.LiquidPhase;
model SaturationFugacitycoefficient2 "fuer phi_sat in Medienmodell"
      // Gmehling und Kolbe, Thermodynamik s.120
  extends BaseFugacityCoefficient2;

input SI.Pressure p_sat[nS];

parameter Integer NoOfEq[nS];

SatFugacityCoeff_TsonopoulosConstants2 satFugacityCoeff_TsonopoulosConstants2(nS=nS, NoOfEq=NoOfEq, Tcrit=Tcrit,mu=mu,pcrit=pcrit);

parameter Real omega[nS];
parameter SI.Temperature Tcrit[nS];
parameter SI.Pressure pcrit[nS];
input Real mu[nS];

protected
Real Tr[nS];
Real f0[nS];
Real f1[nS];
Real f2[nS];
Real f3[nS];
Real Summe[nS];
Real B[nS] "in m³/mol";// zweite Virialkoeffizient des Reinstoffs
Real R = Modelica.Constants.R*1000;
Real pcrit_bar[nS];
//Tsonopoulos Parameter für 2ten Virialkoeff.
Real a[nS]= satFugacityCoeff_TsonopoulosConstants2.a;
Real b[nS]= satFugacityCoeff_TsonopoulosConstants2.b;

Real help1[nS];
Real help2[nS];

equation
  for i in 1:nS loop
  //  omega[i] = Medium.omega[reorgLiq[i]];
    end for;
  for i in 1:nS loop
    pcrit_bar[i] = pcrit[i]/1e5;
    end for;

  for k in 1:nS loop

  //Berechnung des zweiten Virialkoeffizienten
  Tr[k] = T/max(1e-6,Tcrit[k]);
  f0[k] = 0.1445-0.330/Tr[k]-0.1385/Tr[k]^2-0.0121/Tr[k]^3-0.000607/Tr[k]^8;
  f1[k] = 0.0637+0.331/Tr[k]^2-0.423/Tr[k]^3-0.008/Tr[k]^8;
  f2[k] = 1/Tr[k]^6;
  f3[k] = -1/Tr[k]^8;

  Summe[k] = f0[k] + omega[k]*f1[k]+a[k]*f2[k]+b[k]*f3[k];

  B[k] = Modelica.Constants.R*max(1e-6,Tcrit[k])/max(1e-6,pcrit[k])*Summe[k]
      "in m³/mol";

  // Berechnung des Fugazitätskoeffizienten
help1[k]=(1-exp(B[k]*1e3*p_sat[k]/R/T))/3;
help2[k]=exp(B[k]*1e3*p_sat[k]/R/T);
  phi_sat[k]=1-help1[k];
  end for;

end SaturationFugacitycoefficient2;
