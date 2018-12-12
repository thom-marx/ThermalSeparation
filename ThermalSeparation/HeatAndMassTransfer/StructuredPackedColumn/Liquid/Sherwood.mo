within ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Liquid;
model Sherwood "mass transfer coefficient using Sherwood number"
  //extends ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSLiquid;
extends BaseLiqMT;
  Real Sh[n,n_k] "Sherwood number";
  Real Re[n](stateSelect=StateSelect.default) "Reynolds number";
  Real Sc[n,n_k] "Schmidt number";
equation
  for j in 1:n loop
    for i in 1:n_k loop
    Sc[j,i] = max(0,eta[j]/rho[j]./D[j,i]);
    end for;
    Re[j] = max(0,abs(w_sup[j])*geometry.d_char*rho[j]/eta[j]);
    Sh[j,:] = 19*Re[j]^0.5*Sc[j,:].^(1/3);
    k[j,:] =  Sh[j,:].*D[j,:]/geometry.d_char;
  end for;

end Sherwood;
