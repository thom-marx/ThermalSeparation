within ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Vapour;
model RochaStructured "structured packing - Rocha et al."
  //Literature:
  //Rocha et al.:
  //Distillation Columns Containing Structured Packings: A comprehensive Model for Their Performance. 2. Mass Transfer Model
  //Ind. Eng. Chem. Res. 1996, 35, 1660-1667
 // extends ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSVapour;
extends BaseVapMT;
  //  Real Sh[n,a] "Sherwood number";
  Real Re[n](each stateSelect=StateSelect.default) "Reynolds number";
  Real Sc[n,n_k] "Schmidt number";
     parameter Real omega_k=0.05
    "large value if change between constant and variable shall be steep";
   parameter Real omega_time = 200 "Wendepunkt der tanh-Funktion";
   Real omega= 1;// 0.5*tanh(omega_k*(time-omega_time))+0.5;
   parameter Real k_const = 0.06 "constant value to start with";

protected
       SI.Velocity w_eff_l[n] "effective liquid velocity";
     SI.Velocity w_eff_v[n] "effective vapour velocity";
equation
  for j in 1:n loop
    w_eff_l[j] = w_sup_l[j]/(geometry.eps*eps_liq[j]*Modelica.Math.sin(geometry.theta/180*Modelica.Constants.pi));
    w_eff_v[j] = w_sup_v[j]/(geometry.eps*(1-eps_liq[j])*Modelica.Math.sin(geometry.theta/180*Modelica.Constants.pi));
    Sc[j,:] = eta[j]/rho[j]./D[j,:];
    Re[j] = (w_eff_l[j]+w_eff_v[j])*rho[j]*geometry.S/eta[j];
    k[j,:] = fill(k_const,n_k)*(1-omega)+ omega*D[j,:].*Sc[j,:].^(1/3)/geometry.S * 0.054  *Re[j]^(0.8);
    // k_2[j,1] = 0.2*D[j,1]*Sc[j,1].^(1/3)/geometry.S * 0.054  *Re[j]^(0.8);
    // k_2[j,2] =0.2* D[j,2]*Sc[j,2].^(1/3)/geometry.S * 0.054  *Re[j]^(0.8);
    // k_2[j,3] =0.2* D[j,3]*Sc[j,3].^(1/3)/geometry.S * 0.054  *Re[j]^(0.8);
    //     k_2[j,4] = D[j,4]*Sc[j,4].^(1/3)/geometry.S * 0.054  *Re[j]^(0.8);
    //
    //         k_2[j,5] = D[j,5]*Sc[j,5].^(1/3)/geometry.S * 0.054  *Re[j]^(0.8);
    //             k_2[j,6] = D[j,6]*Sc[j,6].^(1/3)/geometry.S * 0.054  *Re[j]^(0.8);
  end for;

  annotation (Documentation(info="<html>
<p>The binary mass transfer coefficients are calculated according to [1]. The correlation is suitable for structured packings.</p>
<p><br/>[1] Rocha, Bravo, Fair: Distillation Columns Containing Structured Packings: A comprehensive model for their performance. 2. Mass Transfer Model, Ind. Eng. Chem. Res., 1996, <i>35</i>, 1660-1667</p>
</html>", revisions="<html>
<pre>Documentation&nbsp;last&nbsp;revised:&nbsp;18.7.2011</pre>
</html>"));
end RochaStructured;
