within ThermalSeparation.Media.IdealGasMixtures;
partial package PartialMediumPure
  "Partial medium properties (base package of all media packages)"

  import      Modelica.Units.SI;
  extends Modelica.Icons.Package;

  // Constants to be set in Medium
  constant String mediumName = "unusablePartialMedium" "Name of the medium";
  constant String substanceNames[:]={mediumName}
    "Names of the mixture substances. Set substanceNames={mediumName} if only one substance.";
  constant String extraPropertiesNames[:]=fill("", 0)
    "Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused";
  constant Boolean singleState
    "= true, if u and d are not a function of pressure";
  constant Boolean reducedX=true
    "= true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance (see docu for details)";
  constant Boolean fixedX=false
    "= true if medium contains the equation X = reference_X";
  constant AbsolutePressure reference_p=101325
    "Reference pressure of Medium: default 1 atmosphere";
  constant Temperature reference_T=298.15
    "Reference temperature of Medium: default 25 deg Celsius";
  constant MassFraction reference_X[nX]= if nX == 0 then fill(0,nX) else fill(1/nX, nX)
    "Default mass fractions of medium";
  constant AbsolutePressure p_default=101325
    "Default value for pressure of medium (for initialization)";
  constant Temperature T_default = Modelica.Units.Conversions.from_degC(20)
    "Default value for temperature of medium (for initialization)";
  constant SpecificEnthalpy h_default = specificEnthalpy_pTX(p_default, T_default, X_default)
    "Default value for specific enthalpy of medium (for initialization)";
  constant MassFraction X_default[nX]=reference_X
    "Default value for mass fractions of medium (for initialization)";

  final constant Integer nS=size(substanceNames, 1) "Number of substances"  annotation(Evaluate=true);
  constant Integer nX=if nS == 1 then 0 else nS
    "Number of mass fractions (= 0, if only one substance)" annotation(Evaluate=true);
  constant Integer nXi=if fixedX then 0 else if reducedX then nS - 1 else nX
    "Number of structurally independent mass fractions (see docu for details)"
   annotation(Evaluate=true);

  final constant Integer nC=size(extraPropertiesNames, 1)
    "Number of extra (outside of standard mass-balance) transported properties"
   annotation(Evaluate=true);

  replaceable record FluidConstants
    "critical, triple, molecular and other standard data of fluid"
    extends Modelica.Icons.Record;
    String iupacName "complete IUPAC name (or common name, if non-existent)";
    String casRegistryNumber
      "chemical abstracts sequencing number (if it exists)";
    String chemicalFormula
      "Chemical formula, (brutto, nomenclature according to Hill";
    String structureFormula "Chemical structure formula";
    MolarMass molarMass "molar mass";
    annotation(Documentation(info="<html></html>"));
  end FluidConstants;

  replaceable record ThermodynamicState
    "Minimal variable set that is available as input argument to every medium function"
    extends Modelica.Icons.Record;
    annotation(Documentation(info="<html></html>"));
  end ThermodynamicState;

  record BasePropertiesRecord
    "Variables contained in every instance of BaseProperties"
    extends Modelica.Icons.Record;
    AbsolutePressure p "Absolute pressure of medium";
    Density d "Density of medium";
    Temperature T "Temperature of medium";
    MassFraction[nX] X(start=reference_X)
      "Mass fractions (= (component mass)/total mass  m_i/m)";
    MassFraction[nXi] Xi(start=reference_X[1:nXi])
      "Structurally independent mass fractions" annotation (Hide=true);
    // SpecificEnthalpy h "Specific enthalpy of medium";
    // SpecificInternalEnergy u "Specific internal energy of medium";
    SpecificHeatCapacity R "Gas constant (of mixture if applicable)";
    MolarMass MM "Molar mass (of mixture or single fluid)";
    annotation(Documentation(info="<html></html>"));
  end BasePropertiesRecord;

  replaceable partial model BaseProperties
    "Base properties (p, d, T, h, u, R, MM and, if applicable, X) of a medium"
    extends BasePropertiesRecord;
    ThermodynamicState state
      "thermodynamic state variables for optional functions";
    parameter Boolean preferredMediumStates=false
      "= true if StateSelect.prefer shall be used for the independent property variables of the medium"
      annotation (Hide=true, Evaluate=true, Dialog(tab="Advanced"));
    Modelica.Units.NonSI.Temperature_degC T_degC=
        Modelica.Units.Conversions.to_degC(T)
      "Temperature of medium in [degC]";
    Modelica.Units.NonSI.Pressure_bar p_bar=
     Modelica.Units.Conversions.to_bar(p)
      "Absolute pressure of medium in [bar]";
    parameter Boolean standardOrderComponents = true
      "if true, last element in components is computed from 1-sum(Xi)";
  equation
    if standardOrderComponents then
      Xi = X[1:nXi];
      if nX > 1 then
        if fixedX then
          X = reference_X;
        end if;
        if reducedX and not fixedX then
          X[nX] = 1 - sum(Xi);
        end if;
        for i in 1:nX loop
          assert(X[i] >= -1.e-5 and X[i] <= 1 + 1.e-5, "Mass fraction X[" +
                 String(i) + "] = " + String(X[i]) + "of substance "
                 + substanceNames[i] + "\nof medium " + mediumName + " is not in the range 0..1");
        end for;
      end if;
    end if;

    assert(p >= 0.0, "Pressure (= " + String(p) + " Pa) of medium \"" +
      mediumName + "\" is negative\n(Temperature = " + String(T) + " K)");

    annotation (structurallyIncomplete,
                Documentation(info="<html>
<p>
Model <b>BaseProperties</b> is a model within package <b>PartialMedium</b>
and contains the <b>declarations</b> of the minimum number of
variables that every medium model is supposed to support.
A specific medium inherits from model <b>BaseProperties</b> and provides
the equations for the basic properties. Note, that in package
PartialMedium the following constants are defined:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td valign=\"top\"><b>Type</b></td>
      <td valign=\"top\"><b>Name</b></td>
      <td valign=\"top\"><b>Description</b></td></tr>
  <tr><td valign=\"top\">String</td><td valign=\"top\">mediumName</td>
      <td valign=\"top\">Unique name of the medium (used to check whether two media in a model
          are the same)</td></tr>
  <tr><td valign=\"top\">String</td><td valign=\"top\">substanceNames</td>
      <td valign=\"top\">Names of the mixture substances that are treated
          as independent.
          If medium consists of a single substance, set substanceNames=fill(\"\",0).
          If medium consists of n substances, provide either n-1 or n
          substance names, depending whether mass fractions
          PartialMedium.BaseProperties.X shall have
          dimension PartialMedium.nX = n-1 or PartialMedium.nX = n</td></tr>
  <tr><td valign=\"top\">Boolean</td><td valign=\"top\">incompressible</td>
      <td valign=\"top\">= true, if density is constant; otherwise set it to false</td></tr>
</table>
<p>
In every medium <b>3+nX equations</b> have to be defined that
provide relations between the following <b>5+nX variables</b>, declared
in model BaseProperties, where nX is the number of independent
mass fractions defined in package PartialMedium:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td valign=\"top\"><b>Variable</b></td>
      <td valign=\"top\"><b>Unit</b></td>
      <td valign=\"top\"><b>Description</b></td></tr>
  <tr><td valign=\"top\">T</td>
      <td valign=\"top\">K</td>
      <td valign=\"top\">temperature</td></tr>
  <tr><td valign=\"top\">p</td>
      <td valign=\"top\">Pa</td>
      <td valign=\"top\">absolute pressure</td></tr>
  <tr><td valign=\"top\">d</td>
      <td valign=\"top\">kg/m^3</td>
      <td valign=\"top\">density</td></tr>
  <tr><td valign=\"top\">h</td>
      <td valign=\"top\">J/kg</td>
      <td valign=\"top\">specific enthalpy</td></tr>
  <tr><td valign=\"top\">u</td>
      <td valign=\"top\">J/kg</td>
      <td valign=\"top\">specific internal energy</td></tr>
  <tr><td valign=\"top\">X[nX]</td>
      <td valign=\"top\">kg/kg</td>
      <td valign=\"top\">independent mass fractions m_i/m</td></tr>
</table>
<p>
In some components, such as \"Ambient\", explicit equations for
medium variables are provided as \"boundary conditions\".
For example, the \"Ambient\" component may define a temperature
T_ambient.
</html>"), Icon(graphics={Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid), Text(extent={{-152,164},{152,
                  102}}, textString =    "%name")}));
  end BaseProperties;

  replaceable partial function setState_pTX
    "Return thermodynamic state as function of p, T and composition X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state "thermodynamic state record";
    annotation(Documentation(info="<html></html>"));
  end setState_pTX;

  replaceable partial function setState_phX
    "Return thermodynamic state as function of p, h and composition X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state "thermodynamic state record";
    annotation(Documentation(info="<html></html>"));
  end setState_phX;

  replaceable partial function setState_psX
    "Return thermodynamic state as function of p, s and composition X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state "thermodynamic state record";
    annotation(Documentation(info="<html></html>"));
  end setState_psX;

  replaceable partial function setState_dTX
    "Return thermodynamic state as function of d, T and composition X or Xi"
    extends Modelica.Icons.Function;
    input Density d "density";
    input Temperature T "Temperature";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state "thermodynamic state record";
    annotation(Documentation(info="<html></html>"));
  end setState_dTX;

  replaceable partial function dynamicViscosity "Return dynamic viscosity"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output DynamicViscosity eta "Dynamic viscosity";
    annotation(Documentation(info="<html></html>"));
  end dynamicViscosity;

  replaceable partial function thermalConductivity
    "Return thermal conductivity"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output ThermalConductivity lambda "Thermal conductivity";
    annotation(Documentation(info="<html></html>"));
  end thermalConductivity;

  replaceable function prandtlNumber "Return the Prandtl number"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output PrandtlNumber Pr "Prandtl number";
  algorithm
    Pr := dynamicViscosity(state)*specificHeatCapacityCp(state)/thermalConductivity(
      state);
  end prandtlNumber;

  replaceable partial function pressure "Return pressure"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output AbsolutePressure p "Pressure";
    annotation(Documentation(info="<html></html>"));
  end pressure;

  replaceable partial function temperature "Return temperature"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output Temperature T "Temperature";
  end temperature;

  replaceable partial function density "Return density"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output Density d "Density";
    annotation(Documentation(info="<html></html>"));
  end density;

  replaceable partial function specificEnthalpy "Return specific enthalpy"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificEnthalpy h "Specific enthalpy";
    annotation(Documentation(info="<html></html>"));
  end specificEnthalpy;

  replaceable partial function specificInternalEnergy
    "Return specific internal energy"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificEnergy u "Specific internal energy";
    annotation(Documentation(info="<html></html>"));
  end specificInternalEnergy;

  replaceable partial function specificEntropy "Return specific entropy"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificEntropy s "Specific entropy";
    annotation(Documentation(info="<html></html>"));
  end specificEntropy;

  replaceable partial function specificGibbsEnergy
    "Return specific Gibbs energy"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificEnergy g "Specific Gibbs energy";
    annotation(Documentation(info="<html></html>"));
  end specificGibbsEnergy;

  replaceable partial function specificHelmholtzEnergy
    "Return specific Helmholtz energy"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificEnergy f "Specific Helmholtz energy";
    annotation(Documentation(info="<html></html>"));
  end specificHelmholtzEnergy;

  replaceable partial function specificHeatCapacityCp
    "Return specific heat capacity at constant pressure"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificHeatCapacity cp
      "Specific heat capacity at constant pressure";
    annotation(Documentation(info="<html></html>"));
  end specificHeatCapacityCp;

  function heatCapacity_cp = specificHeatCapacityCp "alias for deprecated name";
  replaceable partial function specificHeatCapacityCv
    "Return specific heat capacity at constant volume"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificHeatCapacity cv "Specific heat capacity at constant volume";
    annotation(Documentation(info="<html></html>"));
  end specificHeatCapacityCv;

  function heatCapacity_cv = specificHeatCapacityCv "alias for deprecated name";
  replaceable partial function isentropicExponent "Return isentropic exponent"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output IsentropicExponent gamma "Isentropic exponent";
    annotation(Documentation(info="<html></html>"));
  end isentropicExponent;

  replaceable partial function isentropicEnthalpy "Return isentropic enthalpy"
    extends Modelica.Icons.Function;
    input AbsolutePressure p_downstream "downstream pressure";
    input ThermodynamicState refState "reference state for entropy";
    output SpecificEnthalpy h_is "Isentropic enthalpy";
    annotation(Documentation(info="<html></html>"));
  end isentropicEnthalpy;

  replaceable partial function velocityOfSound "Return velocity of sound"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output VelocityOfSound a "Velocity of sound";
    annotation(Documentation(info="<html></html>"));
  end velocityOfSound;

  replaceable partial function isobaricExpansionCoefficient
    "Return overall the isobaric expansion coefficient beta"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output IsobaricExpansionCoefficient beta "Isobaric expansion coefficient";
    annotation(Documentation(info="<html></html>"));
  end isobaricExpansionCoefficient;

  function beta = isobaricExpansionCoefficient
    "alias for isobaricExpansionCoefficient for user convenience";
  replaceable partial function isothermalCompressibility
    "Return overall the isothermal compressibility factor"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SI.IsothermalCompressibility kappa "Isothermal compressibility";
    annotation(Documentation(info="<html></html>"));
  end isothermalCompressibility;

  function kappa = isothermalCompressibility
    "alias of isothermalCompressibility for user convenience";
  // explicit derivative functions for finite element models
  replaceable partial function density_derp_h
    "Return density derivative wrt pressure at const specific enthalpy"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output DerDensityByPressure ddph "Density derivative wrt pressure";
    annotation(Documentation(info="<html></html>"));
  end density_derp_h;

  replaceable partial function density_derh_p
    "Return density derivative wrt specific enthalpy at constant pressure"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output DerDensityByEnthalpy ddhp "Density derivative wrt specific enthalpy";
    annotation(Documentation(info="<html></html>"));
  end density_derh_p;

  replaceable partial function density_derp_T
    "Return density derivative wrt pressure at const temperature"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output DerDensityByPressure ddpT "Density derivative wrt pressure";
    annotation(Documentation(info="<html></html>"));
  end density_derp_T;

  replaceable partial function density_derT_p
    "Return density derivative wrt temperature at constant pressure"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output DerDensityByTemperature ddTp "Density derivative wrt temperature";
    annotation(Documentation(info="<html></html>"));
  end density_derT_p;

  replaceable partial function density_derX
    "Return density derivative wrt mass fraction"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output Density[nX] dddX "Derivative of density wrt mass fraction";
    annotation(Documentation(info="<html></html>"));
  end density_derX;

  replaceable partial function molarMass "Return the molar mass of the medium"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output MolarMass MM "Mixture molar mass";
    annotation(Documentation(info="<html></html>"));
  end molarMass;

  replaceable function specificEnthalpy_pTX
    "Return specific enthalpy from p, T, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:]=reference_X "Mass fractions";
    output SpecificEnthalpy h "Specific enthalpy";
  algorithm
    h := specificEnthalpy(setState_pTX(p,T,X));
    annotation(Documentation(info="<html></html>"));
  end specificEnthalpy_pTX;

  replaceable function density_pTX "Return density from p, T, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:] "Mass fractions";
    output Density d "Density";
  algorithm
    d := density(setState_pTX(p,T,X));
    annotation(Documentation(info="<html></html>"));
  end density_pTX;

  replaceable function temperature_phX
    "Return temperature from p, h, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output Temperature T "Temperature";
  algorithm
    T := temperature(setState_phX(p,h,X));
    annotation(Documentation(info="<html></html>"));
  end temperature_phX;

  replaceable function density_phX "Return density from p, h, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output Density d "Density";
  algorithm
    d := density(setState_phX(p,h,X));
    annotation(Documentation(info="<html></html>"));
  end density_phX;

  replaceable function temperature_psX
    "Return temperature from p,s, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output Temperature T "Temperature";
  algorithm
    T := temperature(setState_psX(p,s,X));
    annotation(Documentation(info="<html></html>"));
  end temperature_psX;

  replaceable function density_psX "Return density from p, s, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output Density d "Density";
  algorithm
    d := density(setState_psX(p,s,X));
    annotation(Documentation(info="<html></html>"));
  end density_psX;

  replaceable function specificEnthalpy_psX
    "Return specific enthalpy from p, s, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output SpecificEnthalpy h "Specific enthalpy";
  algorithm
    h := specificEnthalpy(setState_psX(p,s,X));
    annotation(Documentation(info="<html></html>"));
  end specificEnthalpy_psX;

  type AbsolutePressure = SI.AbsolutePressure (
      min=0,
      max=1.e8,
      nominal=1.e5,
      start=1.e5) "Type for absolute pressure with medium specific attributes";
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
  type SpecificInternalEnergy = SpecificEnergy
    "Type for specific internal energy with medium specific attributes";
  type SpecificEnthalpy = SI.SpecificEnthalpy (
      min=-1.0e8,
      max=1.e8,
      nominal=1.e6)
    "Type for specific enthalpy with medium specific attributes";
  type SpecificEntropy = SI.SpecificEntropy (
      min=-1.e6,
      max=1.e6,
      nominal=1.e3) "Type for specific entropy with medium specific attributes";
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
  package Choices "Types, constants to define menu choices"
    package Init
      "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available"

      extends Modelica.Icons.Package;
      constant Integer NoInit=1;
      constant Integer InitialStates=2;
      constant Integer SteadyState=3;
      constant Integer SteadyMass=4;
      type Temp
        "Temporary type with choices for menus (until enumerations are available)"

        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                NoInit "NoInit (no initialization)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                InitialStates "InitialStates (initialize medium states)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                SteadyState "SteadyState (initialize in steady state)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                SteadyMass
              "SteadyMass (initialize density or pressure in steady state)"));
      end Temp;
      annotation (preferedView="text");
    end Init;

    package ReferenceEnthalpy
      "Type, constants and menu choices to define reference enthalpy, as temporary solution until enumerations are available"

      extends Modelica.Icons.Package;
      constant Integer ZeroAt0K=1;
      constant Integer ZeroAt25C=2;
      constant Integer UserDefined=3;
      type Temp
        "Temporary type with choices for menus (until enumerations are available)"

        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.
                ZeroAt0K
              "The enthalpy is 0 at 0 K (default), if the enthalpy of formation is excluded",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.
                ZeroAt25C
              "The enthalpy is 0 at 25 degC, if the enthalpy of formation is excluded",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.
                UserDefined
              "The user-defined reference enthalpy is used at 293.15 K (25 degC)"));

      end Temp;
      annotation (preferedView="text");
    end ReferenceEnthalpy;

    package ReferenceEntropy
      "Type, constants and menu choices to define reference entropy, as temporary solution until enumerations are available"

      extends Modelica.Icons.Package;
      constant Integer ZeroAt0K=1;
      constant Integer ZeroAt0C=2;
      constant Integer UserDefined=3;
      type Temp
        "Temporary type with choices for menus (until enumerations are available)"

        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                ZeroAt0K "The entropy is 0 at 0 K (default)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                ZeroAt0C "The entropy is 0 at 0 degC",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Init.
                UserDefined
              "The user-defined reference entropy is used at 293.15 K (25 degC)"));

      end Temp;
      annotation (preferedView="text");
    end ReferenceEntropy;

    package pd
      "Type, constants and menu choices to define whether p or d are known, as temporary solution until enumerations are available"

      extends Modelica.Icons.Package;
      constant Integer default=1;
      constant Integer p_known=2;
      constant Integer d_known=3;

      type Temp
        "Temporary type with choices for menus (until enumerations are available)"

        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.pd.default
              "default (no boundary condition for p or d)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.pd.p_known
              "p_known (pressure p is known)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.pd.d_known
              "d_known (density d is known)"));
      end Temp;
      annotation (preferedView="text");
    end pd;

    package Th
      "Type, constants and menu choices to define whether T or h are known, as temporary solution until enumerations are available"

      extends Modelica.Icons.Package;
      constant Integer default=1;
      constant Integer T_known=2;
      constant Integer h_known=3;

      type Temp
        "Temporary type with choices for menus (until enumerations are available)"

        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Th.default
              "default (no boundary condition for T or h)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Th.T_known
              "T_known (temperature T is known)",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Th.h_known
              "h_known (specific enthalpy h is known)"));
      end Temp;
      annotation (preferedView="text");
    end Th;

    package Explicit
      "Type, constants and menu choices to define the explicitly given state variable inputs"

      constant Integer dT_explicit=0 "explicit in density and temperature";
      constant Integer ph_explicit=1
        "explicit in pressure and specific enthalpy";
      constant Integer ps_explicit=2
        "explicit in pressure and specific entropy";
      constant Integer pT_explicit=3 "explicit in pressure and temperature";

      type Temp
        "Temporary type with choices for menus (until enumerations are available)"
        extends Integer(min=0,max=3);
        annotation (Evaluate=true, choices(
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Explicit.dT_explicit
              "explicit in d and T",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Explicit.ph_explicit
              "explicit in p and h",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Explicit.ps_explicit
              "explicit in p and s",
            choice=Modelica.Media.Interfaces.PartialMedium.Choices.Explicit.pT_explicit
              "explicit in p and s"));
      end Temp;
    end Explicit;
  end Choices;
  annotation (Documentation(info="<html>
<p>
<b>PartialMedium</b> is a package and contains all <b>declarations</b> for
a medium. This means that constants, models, and functions
are defined that every medium is supposed to support
(some of them are optional). A medium package
inherits from <b>PartialMedium</b> and provides the
equations for the medium. The details of this package
are described in  
<a href=\"Modelica://Modelica.Media.UsersGuide\">Modelica.Media.UsersGuide</a>.
</p>
</html>
",
 revisions="<html>
  
</html>"));
end PartialMediumPure;
