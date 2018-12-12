within ThermalSeparation.Media.WaterBasedLiquid.N2_CO2_H2O;
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

for i in 1:n loop
  if i==3 then
  h[i] = MM_substances[2]*r_water;
  else
  h[i] = 0;
  end if;
  end for;
  annotation (Icon(graphics={
                        Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
end EvaporationEnthalpy;
