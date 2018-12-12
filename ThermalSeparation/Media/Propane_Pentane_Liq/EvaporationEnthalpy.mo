within ThermalSeparation.Media.Propane_Pentane_Liq;
model EvaporationEnthalpy
/*** enthalpy ***/
input SI.Pressure p;
input SI.Temperature T;
//input SI.MoleFraction x[n];
ThermalSeparation.Units.MolarEnthalpy h[n];

protected
SI.SpecificEnthalpy r_water = -2462.5 * T +3177.8e3;
equation
      /*** enthalpy ***/

  h[1] = MMX[1]*426e3;
  h[2] = MMX[2]*358e3;

  annotation (Icon(graphics={
                        Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
end EvaporationEnthalpy;
