within ThermalSeparation.HeatAndMassTransfer.PackedColumn.Liquid;
model RochaStructured "structured packing - Rocha et al."
 // extends ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSLiquid;
 extends BaseLiqMT;
    //Literature:
  //Rocha et al.:
  //Distillation Columns Containing Structured Packings: A comprehensive Model for Their Performance. 2. Mass Transfer Model
  //Ind. Eng. Chem. Res. 1996, 35, 1660-1667

         SI.Velocity w_eff[n] "effective liquid velocity";
   parameter Real omega_k=0.05
    "large value if change between constant and variable shall be steep";
   parameter Real omega_time = 200 "Wendepunkt der tanh-Funktion";
   Real omega = 1;//0.5*tanh(omega_k*(time-omega_time))+0.5;
   parameter Real k_const = 6e-5 "constant value to start with";
equation
  for j in 1:n loop

    w_eff[j] = w_sup[j]/(geometry.eps*eps_liq[j]*Modelica.Math.sin(geometry.theta/180*Modelica.Constants.pi));
    k[j,:] = fill(k_const,n_k)*(1-omega) + omega *2*(D[j,:]*0.9*w_eff[j]/(Modelica.Constants.pi*geometry.S)).^(0.5);
  end for;

  annotation (Documentation(info="<html>
<p>The binary mass transfer coefficients are calculated according to [1]. The correlation is suitable for structured packings.</p>
<p><br/>[1] Rocha, Bravo, Fair: Distillation Columns Containing Structured Packings: A comprehensive model for their performance. 2. Mass Transfer Model, Ind. Eng. Chem. Res., 1996, <i>35</i>, 1660-1667</p>
</html>", revisions="<html>
<pre>Documentation&nbsp;last&nbsp;revised:&nbsp;18.7.2011</pre>
</html>"));
end RochaStructured;
