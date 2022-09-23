within ThermalSeparation.Media.H2O_CO2_MEA_Liq;
function DensityWater "Return density from p and T"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.AbsolutePressure p "Pressure";
    input Modelica.Units.SI.Temperature T "Temperature";

   output Modelica.Units.SI.Density d "Density";
algorithm
    d :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.rho_pT(
      p, T);
end DensityWater;
