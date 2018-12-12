within ThermalSeparation.Interfaces;
connector LiquidPortOut
  extends Icons.Icons.LiquidPortOut;
 replaceable package Medium = Media.BaseMediumLiquid;
stream SI.MoleFraction x_outflow[Medium.nSubstance];
stream ThermalSeparation.Units.MolarEnthalpy h_outflow;
  flow SI.MolarFlowRate Ndot(max=0, min=-1e12);

  SI.Pressure p;
annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}),
                 graphics={Ellipse(
          extent={{-64,66},{-64,66}},
          lineColor={153,217,234},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          lineThickness=1)}));
end LiquidPortOut;
