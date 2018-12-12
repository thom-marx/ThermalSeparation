within ThermalSeparation.Utilities;
model RefluxConverter "converts r to r/(1+r)"

  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-124,-12},{-76,36}})));
  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{98,-14},{144,32}})));
equation
  y=u/(1+u);

  annotation (Icon(graphics={
        Rectangle(extent={{-100,100},{100,-100}}, lineColor={0,0,255}),
        Text(
          extent={{34,6},{-42,66}},
          lineColor={0,0,255},
          textString="r"),
        Text(
          extent={{38,-58},{-38,2}},
          lineColor={0,0,255},
          textString="r+1"),
        Rectangle(
          extent={{-62,8},{62,4}},
          lineColor={0,0,255},
          lineThickness=1,
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid)}));
end RefluxConverter;
