within ThermalSeparation.Components.SourcesSinks;
model AmbientHeatSink
   extends Icons.Color.AmbientHeatSink;
  ThermalSeparation.Interfaces.HeatPort heatPort1
    annotation (Placement(transformation(extent={{-146,-30},{-106,10}},
                                                                     rotation=0),
        iconTransformation(extent={{-126,-10},{-106,10}})));
  parameter SI.Temperature T_amb = 298 "ambient temperature";
equation
  heatPort1.T = T_amb;
  annotation (Icon(graphics), Diagram(coordinateSystem(preserveAspectRatio=false,
                   extent={{-100,-100},{100,100}}), graphics));
end AmbientHeatSink;
