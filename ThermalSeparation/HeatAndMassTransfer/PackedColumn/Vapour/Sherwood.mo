within ThermalSeparation.HeatAndMassTransfer.PackedColumn.Vapour;
model Sherwood "mass transfer coefficient using Sherwood number"
  //extends ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSVapour;
    extends BaseVapMT;
    Real Sh[n,n_k] "Sherwood number";
  Real Re[n](stateSelect=StateSelect.default) "Reynolds number";
  Real Sc[n,n_k] "Schmidt number";
equation
  for j in 1:n loop
    Sc[j,:] = eta[j]/rho[j]./D[j,:];
    Re[j] = max(0,abs(w_sup_v[j])*geometry.d_char*rho[j]/eta[j]);
    Sh[j,:] = 0.9*geometry.eps^(1/3)*Re[j]^(2/3)*Sc[j,:].^(1/3);
    k[j,:] = Sh[j,:].*D[j,:]/geometry.d_char;
  end for;

end Sherwood;
