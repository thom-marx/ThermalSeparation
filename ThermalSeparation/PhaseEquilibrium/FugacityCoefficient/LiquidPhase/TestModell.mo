within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.LiquidPhase;
model TestModell
parameter Integer n = 1;
parameter Integer nS = 2;

// parameter SI.Temperature T[n] = {343.15};
// parameter SI.Pressure p[n] = {100000};
// parameter SI.Pressure p_sat[n,nS] = {{72300,31090}};
// parameter Real omega[nS] = MediumVapour.omega;
// parameter SI.Temperature Tcrit[nS] = MediumVapour.Tcrit;
// parameter SI.Pressure pcrit[nS] = MediumVapour.pcrit;
// // Wert muss im Medium package wie Tcrit und pcrit ergänzt werden (ThermalSeparation.Media.ModelicaMedia.IdealGases.Common.MixtureGasNasa)
// parameter Integer NoOfEq[nS] = {1,6};

parameter SI.Temperature T[n] = fill(343.15,n);//{343.15};
parameter SI.Pressure p[n] = fill(100000,n);//{100000};
parameter SI.Pressure p_sat[n,nS] = {{72300,31090}};
parameter Real omega[nS] = {Medium.omega[reorgLiq[i]] for i in 1:nS};
parameter SI.Temperature Tcrit[nS] = {Medium.Tcrit[reorgLiq[i]] for i in 1:nS};
parameter SI.Pressure pcrit[nS] = {Medium.pcrit[reorgLiq[i]] for i in 1:nS};
// Wert muss im Medium package wie Tcrit und pcrit ergänzt werden (ThermalSeparation.Media.ModelicaMedia.IdealGases.Common.MixtureGasNasa)
parameter Integer NoOfEq[nS] = {4,6};
parameter Integer reorgLiq[nS] = {1,2};
parameter Boolean inertLiquid[nS]=fill(false,2);

replaceable package Medium = 
      ThermalSeparation.Media.ModelicaMedia.IdealGases.MixtureGases.C2H5OH_H2O 
  constrainedby ThermalSeparation.Media.BaseMediumLiquid;

// SaturationFugacitycoefficient test1(redeclare replaceable package MediumVapour
//       =
//       MediumVapour, n=n, nS=nS, T=T, p=p, p_sat=p_sat, NoOfEq=NoOfEq);

      SaturationFugacitycoefficient test1(redeclare replaceable package Medium
      = 
      Medium, n=n, nS=nS, T=T, p=p, p_sat=p_sat, NoOfEq=NoOfEq);

end TestModell;
