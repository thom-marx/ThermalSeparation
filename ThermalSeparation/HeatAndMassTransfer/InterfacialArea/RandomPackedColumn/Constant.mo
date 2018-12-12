within ThermalSeparation.HeatAndMassTransfer.InterfacialArea.RandomPackedColumn;
model Constant "constant volumetric area"

  extends BasePacked;
  parameter Units.VolumetricArea a_const = geometry.a;
equation
  for j in 1:n loop
 frac[j] = a_const/geometry.a;
  a[j]=a_const;
  end for;
  annotation (Documentation(info="<html>
<p>In this class a constant value for the interfacial area is set.</p>
</html>"));
end Constant;
