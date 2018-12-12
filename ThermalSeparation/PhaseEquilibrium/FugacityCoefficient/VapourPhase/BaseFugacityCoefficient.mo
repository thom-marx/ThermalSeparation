within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase;
partial model BaseFugacityCoefficient "base model for fugacity coefficients"
  parameter Integer n=1 annotation(Dialog(enable=false));
  parameter Integer nS=2 
                       annotation(Dialog(enable=false));
replaceable package MediumVapour = 
      ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2_CO2_SO2_HCl_HF 
                                                            constrainedby
    Media.BaseMediumVapour;

  input SI.Pressure p[n];
  input SI.Temperature T[n];
  input SI.MoleFraction y[n,nS];
  input SI.MolarVolume v[n];
  output Real phi[n,nS];//(start=fill(1,n,nSV));

protected
  Real phi_aux[n,nS];
equation
  /*** phi shall only be calculated for substances which exist in both phases, 
  since for the other substances this information is not needed. 
  The inert substances have an arbitrary value (which is never uses) ***/
//   for i in 1:nS loop
//     phi[:,reorgVap[i]] = phi_aux[:,i];
//     end for;
//
//   for i in 1:nSV loop
//     if inertVapour[i] then
//       phi[:,i] = zeros(n);
//     else
//
//     end if;
//     end for;
phi=phi_aux;
end BaseFugacityCoefficient;
