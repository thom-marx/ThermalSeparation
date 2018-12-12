within ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Vapour;
model MULTIPAK "for MULTIPAK catalytic packing"
  //correlation from Gorak, Hoffmann: Catalytic Distillation in Structured Packings: Methyl Acetate Synthesis, AIChE Journal, 2001, 1067 ff.
  //extends ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSVapour;
    extends BaseVapMT;
    Real Sh[n,n_k] "Sherwood number";
  Real Re[n](stateSelect=StateSelect.default) "Reynolds number";
  Real Sc[n,n_k] "Schmidt number";
protected
  SI.Velocity w_eff_v[n];
  SI.Velocity w_eff_l[n];
   Real omega = 0.5*tanh(0.08*(time-200))+0.5;
equation
  for j in 1:n loop
    w_eff_v[j] = w_sup_v[j]/(geometry.eps*(1-eps_liq[j])*Modelica.Math.sin(geometry.theta));
    w_eff_l[j] = w_sup_l[j]/(geometry.eps*eps_liq[j]*Modelica.Math.sin(geometry.theta));
    Sc[j,:] = eta[j]/rho[j]./D[j,:];
    //Re[j] = max(0,abs(w_sup_v[j])*geometry.S*rho[j]/eta[j]);
     Re[j] = max(0,abs(w_eff_v[j]+w_eff_l[j])*geometry.S*rho[j]/eta[j]);
    Sh[j,:] = 0.0064*Re[j]^0.96*Sc[j,:].^(1/3);
    k[j,:] = fill(0.01,n_k)*(1-omega) + omega*Sh[j,:].*D[j,:]/geometry.S;
  end for;

end MULTIPAK;
