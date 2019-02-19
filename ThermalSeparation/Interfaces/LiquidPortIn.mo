within ThermalSeparation.Interfaces;
connector LiquidPortIn
  extends Icons.Color.LiquidPortIn;
 replaceable package Medium = Media.BaseMediumLiquid;
 stream SI.MoleFraction x_outflow[Medium.nSubstance];
 stream ThermalSeparation.Units.MolarEnthalpy h_outflow;
 flow SI.MolarFlowRate Ndot(min=0,max=1e12);

 SI.Pressure p;

   annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}),
                    graphics={Ellipse(
          extent={{-36,38},{-36,38}},
          lineColor={153,217,234},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          lineThickness=1)}));
end LiquidPortIn;
