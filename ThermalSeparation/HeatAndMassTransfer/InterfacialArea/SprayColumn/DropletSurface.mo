within ThermalSeparation.HeatAndMassTransfer.InterfacialArea.SprayColumn;
model DropletSurface
  parameter Integer n(min=1);
  input SI.Diameter d_drop[n];
  input Real n_drop[n];
  output SI.Area surface[n];

equation
  for j in 1:n loop
  surface[j] = Modelica.Constants.pi*d_drop[j] * d_drop[j]*n_drop[j];
  end for;
  annotation (Documentation(revisions="<html>
<p>The interfacial area is calculated by multiplying the surface of the droplets (assuming them to be perfect spheres) by the number of droplets in one stage.</p>
</html>"));
end DropletSurface;
