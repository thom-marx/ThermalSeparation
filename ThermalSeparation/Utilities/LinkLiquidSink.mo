within ThermalSeparation.Utilities;
model LinkLiquidSink
 // parameter Integer nS(min=1);
 replaceable package Medium=Media.BaseMediumLiquid;
  // parameter Integer n;
  Interfaces.LiquidPortIn liquidIn(redeclare package Medium=Medium);
 // input Modelica.Units.SI.Height z;
 input SI.Pressure p;
equation
  liquidIn.p=p;
 // liquidIn.z = z;
end LinkLiquidSink;
