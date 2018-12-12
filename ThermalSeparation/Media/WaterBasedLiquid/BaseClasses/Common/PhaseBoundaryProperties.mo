within ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.Common;
record PhaseBoundaryProperties
  "thermodynamic base properties on the phase boundary"
  extends Modelica.Icons.Record;
  SI.Density d "density";
  SI.SpecificEnthalpy h "specific enthalpy";
  SI.SpecificEnergy u "inner energy";
  SI.SpecificEntropy s "specific entropy";
  SI.SpecificHeatCapacity cp "heat capacity at constant pressure";
  SI.SpecificHeatCapacity cv "heat capacity at constant volume";
  DerPressureByTemperature pt "derivative of pressure wrt temperature";
  DerPressureByDensity pd "derivative of pressure wrt density";
end PhaseBoundaryProperties;
