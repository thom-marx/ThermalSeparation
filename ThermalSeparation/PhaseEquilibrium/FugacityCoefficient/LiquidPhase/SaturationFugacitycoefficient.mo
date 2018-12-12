within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.LiquidPhase;
model SaturationFugacitycoefficient
      // Gmehling und Kolbe, Thermodynamik s.120
  extends BaseFugacityCoefficient;
replaceable package Medium = 
      ThermalSeparation.Media.WaterBasedLiquid.N2_H2O_O2 
  constrainedby ThermalSeparation.Media.BaseMediumLiquid;
input SI.Pressure p_sat[n,nS];

parameter Integer NoOfEq[nS] = {Medium.eq_Tsonopoulos[reorgLiq[i]] for i in 1:nS};

SatFugacityCoeff_TsonopoulosConstants satFugacityCoeff_TsonopoulosConstants(reorgLiq=reorgLiq, redeclare
      replaceable package Medium = 
      Medium, nS=nS, NoOfEq=NoOfEq);

parameter Real omega[nS] = {Medium.omega[reorgLiq[i]] for i in 1:nS};
parameter SI.Temperature Tcrit[nS]= {Medium.Tcrit[reorgLiq[i]] for i in 1:nS};
parameter SI.Pressure pcrit[nS] = {Medium.pcrit[reorgLiq[i]] for i in 1:nS};

protected
Real Tr[n,nS];
Real f0[n,nS];
Real f1[n,nS];
Real f2[n,nS];
Real f3[n,nS];
Real Summe[n,nS];
Real B[n,nS] "in m³/mol";// zweite Virialkoeffizient des Reinstoffs
Real R = Modelica.Constants.R*1000;
Real pcrit_bar[nS];
//Tsonopoulos Parameter für 2ten Virialkoeff.
Real a[nS]= satFugacityCoeff_TsonopoulosConstants.a;
Real b[nS]= satFugacityCoeff_TsonopoulosConstants.b;

Real help1[n,nS];
Real help2[n,nS];

equation
  for i in 1:nS loop
  //  omega[i] = Medium.omega[reorgLiq[i]];
    end for;
  for i in 1:nS loop
    pcrit_bar[i] = pcrit[i]/1e5;
    end for;
for m in 1:n loop
  for k in 1:nS loop

  //Berechnung des zweiten Virialkoeffizienten
  Tr[m,k] = T[m]/Tcrit[k];
  f0[m,k] = 0.1445-0.330/Tr[m,k]-0.1385/Tr[m,k]^2-0.0121/Tr[m,k]^3-0.000607/Tr[m,k]^8;
  f1[m,k] = 0.0637+0.331/Tr[m,k]^2-0.423/Tr[m,k]^3-0.008/Tr[m,k]^8;
  f2[m,k] = 1/Tr[m,k]^6;
  f3[m,k] = -1/Tr[m,k]^8;

  Summe[m,k] = f0[m,k] + omega[k]*f1[m,k]+a[k]*f2[m,k]+b[k]*f3[m,k];

  B[m,k] = Modelica.Constants.R*Tcrit[k]/pcrit[k]*Summe[m,k] "in m³/mol";

  // Berechnung des Fugazitätskoeffizienten
help1[m,k]=(1-exp(B[m,k]*1e3*p_sat[m,k]/R/T[m]))/3;
help2[m,k]=exp(B[m,k]*1e3*p_sat[m,k]/R/T[m]);
  phi_sat_aux[m,k]=1-help1[m,k];
  end for;
end for;
end SaturationFugacitycoefficient;
