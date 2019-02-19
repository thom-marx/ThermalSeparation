within ThermalSeparation.Interfaces;
connector GasPortOut
  extends Icons.Color.GasPortOut;
  replaceable package Medium = Media.BaseMediumVapour;
 stream SI.MoleFraction x_outflow[Medium.nSubstance];
 stream ThermalSeparation.Units.MolarEnthalpy h_outflow;
SI.Pressure p;
 flow SI.MolarFlowRate Ndot;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}),
                   graphics={Ellipse(
          extent={{-62,68},{-62,68}},
          lineColor={255,127,39},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          lineThickness=1)}));
end GasPortOut;
