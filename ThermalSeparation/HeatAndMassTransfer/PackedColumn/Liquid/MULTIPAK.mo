within ThermalSeparation.HeatAndMassTransfer.PackedColumn.Liquid;
model MULTIPAK "for MULTIPAK catalytic packing"
  //correlation from :
  //Gorak, Hoffmann: Catalytic Distillatin in Structured Packings: Methyl Acetate Synthesis, AIChE Journal, Vol. 47, No. 5, 2001, p. 1067- ...
extends BaseLiqMT;

protected
  SI.Velocity w_eff_l[n];
equation

  for j in 1:n loop
     w_eff_l[j] = max(1e-3,w_sup[j]/(geometry.eps*eps_liq[j]*Modelica.Math.sin(geometry.theta)));
    for i in 1:n_k loop
    k[j,i] = 2*sqrt(0.9*D[j,i]*w_eff_l[j]/(Modelica.Constants.pi*geometry.S));
    end for;
  end for;

end MULTIPAK;
