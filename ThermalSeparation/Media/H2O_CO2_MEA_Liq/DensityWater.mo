within ThermalSeparation.Media.H2O_CO2_MEA_Liq;
function DensityWater "Return density from p and T"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.AbsolutePressure p "Pressure";
    input Modelica.SIunits.Temperature T "Temperature";

   output Modelica.SIunits.Density d "Density";
algorithm
    d :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.rho_pT(
      p, T);
end DensityWater;
