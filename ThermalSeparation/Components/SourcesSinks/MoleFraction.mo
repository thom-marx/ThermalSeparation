within ThermalSeparation.Components.SourcesSinks;
model MoleFraction "source for mole fractions"
parameter Integer nS(min=1)=1;
parameter Real duration= 10;
parameter Real startTime = 1;
parameter Real startValue[nS];
parameter Real endValue[nS];
  Modelica.Blocks.Interfaces.RealOutput y[nS]
    annotation (Placement(transformation(extent={{70,-8},{108,30}}),
        iconTransformation(extent={{70,-8},{108,30}})));
  Modelica.Blocks.Sources.Ramp[nS] ramp(each duration=duration, each
      startTime =                                                              startTime,offset=startValue,height=endValue-startValue);
equation
  for i in 1:nS loop
    y[i]=ramp[i].y;
  end for;

  annotation (Icon(graphics={Rectangle(
          extent={{-60,80},{80,-60}},
          fillColor={135,135,135},
          fillPattern=FillPattern.Sphere,
          pattern=LinePattern.None,
          lineColor={0,0,0}), Bitmap(extent={{108,12},{104,18}}, fileName="")}));
end MoleFraction;
