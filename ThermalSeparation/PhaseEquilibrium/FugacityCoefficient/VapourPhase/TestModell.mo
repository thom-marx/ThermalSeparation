within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase;
model TestModell
parameter Integer n = 1;
parameter Integer nS = 2;

parameter SI.Temperature T[n] = {144.26};
parameter SI.Pressure p[n] = {100000};
parameter SI.Pressure psat[nS] = {31090,72300};//,7214000};
parameter SI.MolarVolume v[n] = {0.044/1.7};
parameter SI.MoleFraction x[n,nS] = {{0.85,0.1}};//,0.05}};
parameter SI.MoleFraction y[n,nS] = {{0.1,0.5}};//,0.4}};
parameter Real omega[nS] = MediumVapour.omega;
parameter SI.Temperature Tcrit[nS] = MediumVapour.Tcrit;
parameter SI.Pressure pcrit[nS] = MediumVapour.pcrit;

  replaceable package MediumVapour = 
      ThermalSeparation.Media.IdealGasMixtures.N2_H2O;
      //für N2-CH4 Mischung: kij = 0.0267

//  Virial_SecondVirialCoefficient Virial_Coeff(redeclare replaceable package
//       MediumVapour =
//             MediumVapour,n=n, nS=nS, T=T);

 Virial virial(                      redeclare replaceable package MediumVapour
      =    MediumVapour,n=n, nS=2, T = fill(144.26,n), p = fill(2000000,n),x = fill({0.2152, 0.7848},n), y = fill({0.6, 0.4},n), v=fill(0.0004601,n));

   ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase.SoaveRedlichKwong
    SRK(
    redeclare replaceable package MediumVapour = MediumVapour,
    n=n,
    nS=2,
    T=fill(144.26, n),
    p=fill(2000000, n),
    x=fill({0.2152,0.7848}, n),
    y=fill({0.6,0.4}, n),
    v=fill(0.0004601, n));

  ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase.PengRobinson
    PR(
    redeclare replaceable package MediumVapour = MediumVapour,
    n=n,
    nS=2,
    T=fill(144.26, n),
    p=fill(2000000, n),
    x=fill({0.2152,0.7848}, n),
    y=fill({0.6,0.4}, n),
    v=fill(0.0004601, n));

end TestModell;
