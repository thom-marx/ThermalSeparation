within ThermalSeparation.Utilities;
model LinkVapourSink
// parameter Integer nS(min=1);
replaceable package Medium=Media.BaseMediumVapour;
 // parameter Integer n( min=1);
 Interfaces.GasPortIn vapourIn(redeclare package Medium=Medium);
 input SI.Pressure p;
equation
 vapourIn.p = p;
end LinkVapourSink;
