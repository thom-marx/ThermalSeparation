within ThermalSeparation;
package Units "Type Definitions"
  import      Modelica.Units.SI;

  type CoefficentOfMassTransfer = Real (quantity="CoefficentOfMassTransfer", final unit="m/s");
  type F_Factor =Real (quantity="F_Factor", final unit="1");
  type MolarDensity = Real (quantity="MolarDensity", final unit =   "mol/m3");
  type MolarEnergy =         Real (quantity="MolarSpecificEnergy", final unit =   "J/mol");
  type MolarEnthalpy =         ThermalSeparation.Units.MolarEnergy;
  type ReactionRate =  SI.MolarFlowRate;
  type VolumetricArea =Real (quantity="VolumetricArea", final unit="m2/m3");
  type DipoleMoment = Real (
    min=0.0,
    max=2.0,
    unit="debye",
    quantity="ElectricDipoleMoment")
    "dipole moment with non-SI-unit debye as returned by the medium model";
  type ConversionDebyeCm = Real (unit="debye/(C.m)",
    quantity="ConversionFactor")
    "dipole moment with non-SI-unit debye as returned by the medium model";
annotation(preferedView="info");
end Units;
