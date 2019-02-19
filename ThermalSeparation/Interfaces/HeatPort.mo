within ThermalSeparation.Interfaces;
connector HeatPort
  extends Icons.Color.HeatPort;
 //parameter Integer n;
 flow SI.HeatFlowRate Qdot;
 SI.Temperature T;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={Ellipse(
          extent={{-24,28},{-24,28}},
          lineColor={188,51,69},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}));
end HeatPort;
