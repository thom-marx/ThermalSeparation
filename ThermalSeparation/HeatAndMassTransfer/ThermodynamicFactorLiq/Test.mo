within ThermalSeparation.HeatAndMassTransfer.ThermodynamicFactorLiq;
model Test
  parameter Integer nS=5;
  parameter Integer n=3;
  parameter Real x[n,nS]=fill(ones(nS)*1/nS,n);
  //parameter Real R=10;
  parameter Real T[n]={400,300,200};
  //parameter Real V[n,nS]={{10,100},{10,100},{100,10}};
  ThermalSeparation.HeatAndMassTransfer.ThermodynamicFactorLiq.NRTLMulticomponent
    nrtlMulticomponent(
    nS=nS,
    n=n,
    x=x,
    T=T);
  WilsonMulticomponent wilsonMulticomponent(nS=nS, n=n, x=x, T=T);
  ThermalSeparation.HeatAndMassTransfer.ThermodynamicFactorLiq.UNIQUACMulticomponent
    uniquacMulticomponent(
    nS=nS,
    n=n,
    x=x,
    T=T);
  //UniquacBinary UniquacBinary(nS=nS, n=n,x=x,T=T);
  Ideal ideal(nS=nS, n=n,x=x,T=T);
equation

end Test;
