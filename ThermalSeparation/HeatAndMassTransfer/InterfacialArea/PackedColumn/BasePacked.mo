within ThermalSeparation.HeatAndMassTransfer.InterfacialArea.PackedColumn;
partial model BasePacked
  extends Base;
      replaceable record Geometry = 
      ThermalSeparation.Geometry.PackedColumn.Geometry                            constrainedby
    ThermalSeparation.Geometry.PackedColumn.Geometry;
    Geometry geometry(final n=n);
    parameter Real frickel=1 "factor to adjust results to measurment data";

protected
     SI.Velocity w_sup[n] "superficial liquid velocity";
     Real frac[n] "effective interfacial area over geometric interfacial area";
equation
  for j in 1:n loop
  w_sup[j] = max(1e-10,Vdot_l[j]/geometry.A);
  end for;

end BasePacked;
