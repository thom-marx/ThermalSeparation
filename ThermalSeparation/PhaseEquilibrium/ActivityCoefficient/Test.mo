within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient;
model Test
  parameter Integer n=4;
  parameter Integer nS=3;
  Wilson Wilson_eq(nS=nS, x_l={1,0,0}, T=300);
 NRTL NRTL_eq(nS=nS, x_l=fill(0.5,nS), T=200);
  UNIQUAC UNIQUAC_eq( nS=nS, x_l=fill(0.5,nS), T=200);
 Margules Margules_eq( nS=nS, x_l=fill(0.5,nS), T=200);//nochmal gucken... nach dem a
equation

end Test;
