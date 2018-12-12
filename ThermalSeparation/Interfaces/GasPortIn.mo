within ThermalSeparation.Interfaces;
connector GasPortIn
  extends Icons.Icons.GasPortIn;
 replaceable package Medium = Media.BaseMediumVapour;
  stream SI.MoleFraction x_outflow[Medium.nSubstance];
  stream ThermalSeparation.Units.MolarEnthalpy h_outflow;
 flow SI.MolarFlowRate Ndot;
SI.Pressure p;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                       graphics={Ellipse(
          extent={{-20,22},{-20,22}},
          lineColor={255,127,39},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          lineThickness=1)}));
end GasPortIn;
