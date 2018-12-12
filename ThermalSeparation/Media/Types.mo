within ThermalSeparation.Media;
package Types
  type AbsolutePressure = SI.AbsolutePressure (
        min=0,
        max=1.e8,
        nominal=1.e5,
        start=1.e5)
    "Type for absolute pressure with medium specific attributes";
  type Density = SI.Density (
      min=0,
      max=1.e5,
      nominal=1,
      start=1) "Type for density with medium specific attributes";
  type DynamicViscosity = SI.DynamicViscosity (
      min=0,
      max=1.e8,
      nominal=1.e-3,
      start=1.e-3) "Type for dynamic viscosity with medium specific attributes";
  type EnthalpyFlowRate = SI.EnthalpyFlowRate (
      nominal=1000.0,
      min=-1.0e8,
      max=1.e8) "Type for enthalpy flow rate with medium specific attributes";
  type MassFlowRate = SI.MassFlowRate (
      quantity="MassFlowRate." + mediumName,
      min=-1.0e5,
      max=1.e5) "Type for mass flow rate with medium specific attributes";
  type MassFraction = Real (
      quantity="MassFraction",
      final unit="kg/kg",
      min=0,
      max=1,
      nominal=0.1) "Type for mass fraction with medium specific attributes";
  type MoleFraction = Real (
      quantity="MoleFraction",
      final unit="mol/mol",
      min=0,
      max=1,
      nominal=0.1) "Type for mole fraction with medium specific attributes";
  type MolarMass = SI.MolarMass (
      min=0.001,
      max=0.25,
      nominal=0.032) "Type for molar mass with medium specific attributes";
  type MolarVolume = SI.MolarVolume (
      min=1e-6,
      max=1.0e6,
      nominal=1.0) "Type for molar volume with medium specific attributes";
  type IsentropicExponent = SI.RatioOfSpecificHeatCapacities (
      min=1,
      max=500000,
      nominal=1.2,
      start=1.2) "Type for isentropic exponent with medium specific attributes";
  type SpecificEnergy = SI.SpecificEnergy (
      min=-1.0e8,
      max=1.e8,
      nominal=1.e6) "Type for specific energy with medium specific attributes";
  type SpecificInternalEnergy = ThermalSeparation.Media.Types.SpecificEnergy
    "Type for specific internal energy with medium specific attributes";
  type SpecificEntropy = SI.SpecificEntropy (
      min=-1.e6,
      max=1.e6,
      nominal=1.e3) "Type for specific entropy with medium specific attributes";
  type SpecificEnthalpy = SI.SpecificEnthalpy (
      min=-1.0e8,
      max=1.e8,
      nominal=1.e6)
    "Type for specific enthalpy with medium specific attributes";
  type SpecificHeatCapacity = SI.SpecificHeatCapacity (
      min=0,
      max=1.e6,
      nominal=1.e3,
      start=1.e3)
    "Type for specific heat capacity with medium specific attributes";
  type SurfaceTension = SI.SurfaceTension
    "Type for surface tension with medium specific attributes";
  type Temperature = SI.Temperature (
      min=1,
      max=1.e4,
      nominal=300,
      start=300) "Type for temperature with medium specific attributes";
  type ThermalConductivity = SI.ThermalConductivity (
      min=0,
      max=500,
      nominal=1,
      start=1) "Type for thermal conductivity with medium specific attributes";
  type PrandtlNumber = SI.PrandtlNumber (
      min=1e-3,
      max=1e5,
      nominal=1.0) "Type for Prandtl number with medium specific attributes";
  type VelocityOfSound = SI.Velocity (
      min=0,
      max=1.e5,
      nominal=1000,
      start=1000) "Type for velocity of sound with medium specific attributes";
  type ExtraProperty = Real (min=0.0, start=1.0)
    "Type for unspecified, mass-specific property transported by flow";
  type CumulativeExtraProperty = Real (min=0.0, start=1.0)
    "Type for conserved integral of unspecified, mass specific property";
  type ExtraPropertyFlowRate = Real
    "Type for flow rate of unspecified, mass-specific property";
  type IsobaricExpansionCoefficient = Real (
      min=1e-8,
      max=1.0e8,
      unit="1/K")
    "Type for isobaric expansion coefficient with medium specific attributes";
  type DipoleMoment = Real (
      min=0.0,
      max=2.0,
      unit="debye",
      quantity="ElectricDipoleMoment")
    "Type for dipole moment with medium specific attributes";
  type DerDensityByPressure = SI.DerDensityByPressure
    "Type for partial derivative of density with resect to pressure with medium specific attributes";
  type DerDensityByEnthalpy = SI.DerDensityByEnthalpy
    "Type for partial derivative of density with resect to enthalpy with medium specific attributes";
  type DerEnthalpyByPressure = SI.DerEnthalpyByPressure
    "Type for partial derivative of enthalpy with resect to pressure with medium specific attributes";
  type DerDensityByTemperature = SI.DerDensityByTemperature
    "Type for partial derivative of density with resect to temperature with medium specific attributes";
  type FixedPhase = Integer(min=0,max=2)
    "phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g. interactive use";
end Types;
