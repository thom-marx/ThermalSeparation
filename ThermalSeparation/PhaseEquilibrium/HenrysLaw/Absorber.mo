within ThermalSeparation.PhaseEquilibrium.HenrysLaw;
model Absorber "frickel solution for Absorber: CO2, H2O"
  extends BaseHenry;

equation
  for j in 1:n loop
 for i in 1:nSV loop

     He[j,i] = if henry_temp then 
                 exp(21.2220117 - 8123.6037/T[j])*1e5 else 
                 1e4;
 end for;
 end for;
end Absorber;
