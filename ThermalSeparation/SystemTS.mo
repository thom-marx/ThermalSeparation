within ThermalSeparation;
model SystemTS "system model"
 parameter SI.Temperature T_ref = 293.15 "reference temperature";

 annotation (
    defaultComponentName="systemTS",
    defaultComponentPrefixes="inner",
    Diagram(graphics),
    Icon(graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={112,146,190},
          fillColor={112,146,190},
          fillPattern=FillPattern.Solid), Rectangle(
          extent={{-80,80},{80,-80}},
          lineColor={112,146,190},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-56,14},{60,-26}},
          lineColor={0,0,0},
          textString="%name")}));

end SystemTS;
