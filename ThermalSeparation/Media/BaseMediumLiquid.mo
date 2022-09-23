within ThermalSeparation.Media;
package BaseMediumLiquid "base class for liquid medium package"
   constant Integer nSubstance(min=2)=2 "number of components";
   constant SI.Temperature Tcrit[nSubstance];
   constant SI.Pressure pcrit[nSubstance];

   constant SI.MolarVolume Vcrit[nSubstance];
   constant Real omega[nSubstance] "acentric factor";
   constant SI.ElectricDipoleMomentOfMolecule mu[nSubstance];
   constant SI.MolarMass MMX[nSubstance];
   constant Integer eq_Tsonopoulos[nSubstance]
    "no. of equation for equation a and b in Table 4-5 in Properties of Gases and Liquids, 5th ed.";
   constant Boolean henry[nSubstance]
    "true if liquid fugacity of substance is to be calculated with Henry's law";
   constant Boolean has_etaSubstance[nSubstance]
    "true if substance i provides a value for the pure dynamic viscosity at system pressure and temperature";
   constant Integer ic[nSubstance] =  zeros(nSubstance)
    "ionic charge of the components, i.e. H+ = 1, OH- = -1, H2O=0";

record ThermodynamicState "thermodynamic state"
 extends Modelica.Icons.Record;
  ThermalSeparation.Media.Types.SpecificEnthalpy h "specific enthalpy";
  ThermalSeparation.Media.Types.Density d "density";
  ThermalSeparation.Media.Types.Temperature T "temperature";
  ThermalSeparation.Media.Types.AbsolutePressure p "pressure";
  ThermalSeparation.Media.Types.FixedPhase phase(
                   min=0, max=2)= 1;
end ThermodynamicState;

replaceable partial model BaseProperties
    "Base properties (p, d, T, h, u, R, MM and, if applicable, X) of a medium"

  input ThermalSeparation.Media.Types.AbsolutePressure p
      "Absolute pressure of medium";
  input ThermalSeparation.Media.Types.Temperature T "Temperature of medium";
  input SI.MoleFraction x[nSubstance];
  output SI.Concentration c[nSubstance] = x*d/MM;

  parameter SI.Temperature T0 = 293.15 "reference temperature";

  ThermodynamicState state
      "thermodynamic state variables for optional functions";

  ThermodynamicProperties properties;

  output SI.SurfaceTension sigma;
  output SI.DynamicViscosity eta;
  output SI.Pressure p_sat[nSubstance] "saturation pressure of the components";

  input ThermalSeparation.Units.MolarEnthalpy h "Specific enthalpy of medium";
  output SI.MolarInternalEnergy u "Specific internal energy of medium";
  output ThermalSeparation.Media.Types.Density d "Density of medium";

  output ThermalSeparation.Media.Types.MolarMass MM "Molar mass of mixture";
  output SI.MolarVolume v(start=1e-5,stateSelect=StateSelect.never); //unit: m3/mol
  output SI.ThermalConductivity lambda;
  output SI.SpecificHeatCapacity cp;

    //dynamic viscosity for each component in the liquid phase
  //if a pure component does not exist as liquid at system pressure and temperature any value can be given
  //as long as the corresponding entry in the constant eta_pure is false
  output SI.DynamicViscosity eta_comp[nSubstance];

  /*** for feed stage ***/
  ThermalSeparation.Interfaces.MediumConOut mediumConOut(h=h, rho = d, MM=MM);

equation
  state.p = p;
  state.T = T;
  state.d = d;
  state.h = h/0.018;

  properties.T=T;
  properties.eta=eta;
  properties.sigma=sigma;
  properties.rho=d;
  properties.MM = MM;
  properties.v = v;
  properties.x = x;
  properties.lambda = lambda;
  properties.cp = cp;
  properties.u = u;
  properties.p = p;
  properties.eta_comp=eta_comp;
  properties.c = c;
  properties.d=d;
  properties.h=h;

  annotation (Icon(graphics={
                        Rectangle(
            extent={{-102,100},{98,-100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
end BaseProperties;

replaceable partial model ActivityCoefficient
  input SI.Temperature T;
  input SI.MoleFraction x_l[nSubstance];
  output Real gamma[nSubstance] "activity coefficient";

  annotation (Icon(graphics={
                      Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}));
end ActivityCoefficient;

replaceable partial model FugacityCoefficient
    "calculates the fugacity coefficient at the saturation point for each component in the liquid phase"
input SI.Pressure p;
input SI.Temperature T;
input SI.Pressure p_sat[nSubstance];
output Real phi_sat[nSubstance] "fugacity coefficient of the saturated vapour";

  annotation (Icon(graphics={
                      Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}));
end FugacityCoefficient;

replaceable model HenryCoefficient
    "calculates Henry coefficient of each substance in the liquid phase, if applicable, otherwise dummy values are used"

  input SI.MoleFraction x_l[nSubstance];
  input SI.Temperature T;
  output SI.Pressure He[nSubstance];

equation
  /*** dummy value for Henry coefficient assigned, can be used if liquid fugacity is not calculated using the Henry coefficient ***/
for i in 1:nSubstance loop
  if not henry[i] then
  He[i] = 0;
  end if;
end for;
  annotation (Icon(graphics={
                      Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}));
end HenryCoefficient;

replaceable partial model DiffusionCoefficient
    "returns the binary diffusion coefficients of the components in the mixture"

   input SI.Temperature T;
   input SI.Pressure p;
   input SI.MoleFraction x[nSubstance];
   input SI.DynamicViscosity eta[nSubstance]
      "viscosity of pure liquid at system pressure and temperature";

   // Diffusion coefficents in the following form: D12, D13, D14, D23, D24, D34
   // the numbers correspond to the ordering of the substances
   output SI.DiffusionCoefficient D[a](start=fill(1e-5,a)) "Maxwell-Stefan diffusion coefficients";
   output SI.DiffusionCoefficient D_matrix[nSubstance,nSubstance]
      "Maxwell-Stefan diffusion coefficients in matrix form, where D_12 = D_21";

   //for diluted systems:
   parameter Boolean diluted = false
      "system is diluted: if the variable diluted is set to true, the extending class must supply an equation for D_diluted";
   parameter Integer solvent = 1 "index of solvent for diluted systems";
   output SI.DiffusionCoefficient D_diluted[nSubstance]
      "diffusion coefficients of the binary system consisting of the particular species and the solvent";

  protected
  parameter Integer aux[:] = {1,3,6,10,15, 21, 28, 36, 45};
  parameter Integer a = aux[nSubstance-1]
      "number of binary diffusion coefficients depending on the number of substances";
equation
 if not diluted then
   D_diluted = fill(0,nSubstance)
        "dummy values for the variable D_diluted; not to be used";
 end if;

     for i in 1:nSubstance loop
       for m in i:nSubstance loop
         if i==m then
           //dummy values for entries on the diagonal, since they are never used
           //To Do: put self-diffusion coefficients here
           D_matrix[i,m] = 10;
         else
           //entries of vector D are put on the right place in the matrix
           D_matrix[i,m] = D[m-2+i];
           end if;
         end for;
         for m in 1:i-1 loop
           //D_12 = D_21, D_13 = D_31 etc.
           D_matrix[i,m] = D_matrix[m,i];
         end for;
     end for;

  annotation (Icon(graphics={
                      Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}));
end DiffusionCoefficient;

  replaceable model SaturationPressure
    Real[8] o_s "vector of auxiliary variables";
      input ThermalSeparation.Media.Types.Temperature T "Temperature of medium";
      output SI.Pressure p_sat[nSubstance]
      "Sttigungspartialdruck der reinen Komponenten";
    Real Tlim=min(T, 647.096);

  end SaturationPressure;

  replaceable partial model CalcSpecificEnthalpy
  /*** enthalpy ***/
  parameter SI.Temperature T0 = 293.15 "reference temperature";
  input SI.Pressure p;
  input SI.Temperature T;
  input SI.MoleFraction x[nSubstance];
  output ThermalSeparation.Units.MolarEnthalpy h;

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;

  replaceable partial model ThermodynamicFactor
    constant Real R=Modelica.Constants.R;
    input SI.Temperature T;
    input SI.MoleFraction x[nSubstance];
    output Real Gamma[nSubstance,nSubstance];
  equation

    annotation (Icon(graphics={
                        Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end ThermodynamicFactor;

record ThermodynamicProperties "thermodynamic properties"
 extends Modelica.Icons.Record;

  ThermalSeparation.Media.Types.Temperature T "temperature";

  SI.SurfaceTension sigma;
  SI.DynamicViscosity eta;
  SI.DynamicViscosity eta_comp[nSubstance];

  ThermalSeparation.Media.Types.Density rho "Density of medium";

  ThermalSeparation.Media.Types.MolarMass MM
      "Molar mass (of mixture or single fluid)";
  SI.MolarVolume v(start=1e-5,stateSelect=StateSelect.never); //Einheit: m3/mol
  SI.MoleFraction x[nSubstance];
  SI.Density d;
  ThermalSeparation.Units.MolarEnthalpy h;

  SI.ThermalConductivity lambda;
  SI.SpecificHeatCapacity cp;
  SI.MolarInternalEnergy u;
  SI.Pressure p;
  SI.Concentration c[nSubstance];

end ThermodynamicProperties;

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
      nominal=300) "Type for temperature with medium specific attributes";
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
end BaseMediumLiquid;
