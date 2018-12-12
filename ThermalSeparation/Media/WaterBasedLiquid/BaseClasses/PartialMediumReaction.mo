within ThermalSeparation.Media.WaterBasedLiquid.BaseClasses;
partial package PartialMediumReaction
  "Partial medium properties (base package of all media packages)"
  extends ThermalSeparation.Media.BaseMediumLiquidReaction;
  import SI = Modelica.SIunits;
  extends Modelica.Icons.Library;

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
  constant Temperature T_default = Modelica.SIunits.Conversions.from_degC(20)
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

  constant Boolean ph_explicit= false
    "true if explicit in pressure and specific enthalpy";
  constant Boolean dT_explicit = false
    "true if explicit in density and temperature";
  constant Boolean pT_explicit= true
    "true if explicit in pressure and temperature";

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
        Temperature criticalTemperature "critical temperature";
      AbsolutePressure criticalPressure "critical pressure";
      MolarVolume criticalMolarVolume "critical molar Volume";
      Real acentricFactor "Pitzer acentric factor";
      Temperature triplePointTemperature "triple point temperature";
      AbsolutePressure triplePointPressure "triple point pressure";
      Temperature meltingPoint "melting point at 101325 Pa";
      Temperature normalBoilingPoint "normal boiling point (at 101325 Pa)";
      DipoleMoment dipoleMoment
      "dipole moment of molecule in Debye (1 debye = 3.33564e10-30 C.m)";
      Boolean hasIdealGasHeatCapacity=false
      "true if ideal gas heat capacity is available";
      Boolean hasCriticalData=false "true if critical data are known";
      Boolean hasDipoleMoment=false "true if a dipole moment known";
      Boolean hasFundamentalEquation=false "true if a fundamental equation";
      Boolean hasLiquidHeatCapacity=false
      "true if liquid heat capacity is available";
      Boolean hasSolidHeatCapacity=false
      "true if solid heat capacity is available";
      Boolean hasAccurateViscosityData=false
      "true if accurate data for a viscosity function is available";
      Boolean hasAccurateConductivityData=false
      "true if accurate data for thermal conductivity is available";
      Boolean hasVapourPressureCurve=false
      "true if vapour pressure data, e.g. Antoine coefficents are known";
      Boolean hasAcentricFactor=false
      "true if Pitzer accentric factor is known";
      SpecificEnthalpy HCRIT0=0.0
      "Critical specific enthalpy of the fundamental equation";
      SpecificEntropy SCRIT0=0.0
      "Critical specific entropy of the fundamental equation";
      SpecificEnthalpy deltah=0.0
      "Difference between specific enthalpy model (h_m) and f.eq. (h_f) (h_m - h_f)";
      SpecificEntropy deltas=0.0
      "Difference between specific enthalpy model (s_m) and f.eq. (s_f) (s_m - s_f)";

    annotation(Documentation(info="<html></html>"));
  end FluidConstants;

constant FluidConstants[nS] fluidConstants "constant data for the fluid";

  replaceable record SaturationProperties
    "Saturation properties of two phase medium"
    extends Modelica.Icons.Record;
    AbsolutePressure psat "saturation pressure";
    Temperature Tsat "saturation temperature";
    annotation(Documentation(info="<html></html>"));
  end SaturationProperties;

  redeclare replaceable model extends BaseProperties
    "Base properties (p, d, T, h, u, R, MM and, if applicable, X) of a medium"
    MassFraction[nX] X(start=reference_X)
      "Mass fractions (= (component mass)/total mass  m_i/m)";
    MassFraction[nXi] Xi(start=reference_X[1:nXi])
      "Structurally independent mass fractions" annotation (Hide=true);
    SpecificHeatCapacity R "Gas constant (of mixture if applicable)";

    parameter Boolean preferredMediumStates=false
      "= true if StateSelect.prefer shall be used for the independent property variables of the medium"
      annotation (Hide=true, Evaluate=true, Dialog(tab="Advanced"));
    SI.Conversions.NonSIunits.Temperature_degC T_degC=
        Modelica.SIunits.Conversions.to_degC(T)
      "Temperature of medium in [degC]";
    SI.Conversions.NonSIunits.Pressure_bar p_bar=
     Modelica.SIunits.Conversions.to_bar(p)
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

 replaceable function setState_pT "Return thermodynamic state from p and T"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
      input ThermalSeparation.Media.Types.FixedPhase phase=
                             0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state "thermodynamic state record";
 algorithm
    state := setState_pTX(p,T,fill(0,0),phase);
 end setState_pT;

  replaceable function setState_ph "Return thermodynamic state from p and h"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
      input ThermalSeparation.Media.Types.FixedPhase phase=
                             0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state "thermodynamic state record";
  algorithm
    state := setState_phX(p,h,fill(0, 0),phase);
  end setState_ph;

  replaceable function setState_ps "Return thermodynamic state from p and s"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
      input ThermalSeparation.Media.Types.FixedPhase phase=
                             0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state "thermodynamic state record";
  algorithm
    state := setState_psX(p,s,fill(0,0),phase);
  end setState_ps;

  replaceable function setState_dT "Return thermodynamic state from d and T"
    extends Modelica.Icons.Function;
    input Density d "density";
    input Temperature T "Temperature";
      input ThermalSeparation.Media.Types.FixedPhase phase=
                             0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state "thermodynamic state record";
  algorithm
    state := setState_dTX(d,T,fill(0,0),phase);
  end setState_dT;

  function density_ph "Return density from p and h"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
   input ThermalSeparation.Media.Types.FixedPhase phase=
                          0 "2 for two-phase, 1 for one-phase, 0 if not known";

    output Density d "Density";
  algorithm
   // d := density_phX(p, h, fill(0,0), phase);
     d :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.rho_ph(
        p,
        h,
        phase);
  end density_ph;

  function temperature_ph "Return temperature from p and h"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
      input ThermalSeparation.Media.Types.FixedPhase phase=
                             0
      "2 for two-phase, 1 for one-phase, 0 if not known";

    output Temperature T "Temperature";
  algorithm
    T :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.T_ph(
        p,
        h,
        phase);
  end temperature_ph;

  function pressure_dT "Return pressure from d and T"
    extends Modelica.Icons.Function;
    input Density d "Density";
    input Temperature T "Temperature";
      input ThermalSeparation.Media.Types.FixedPhase phase=
                             0
      "2 for two-phase, 1 for one-phase, 0 if not known";

    output AbsolutePressure p "Pressure";
  algorithm
     p :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.p_dT(
        d,
        T,
        phase);
  end pressure_dT;

  function specificEnthalpy_dT "Return specific enthalpy from d and T"
    extends Modelica.Icons.Function;
    input Density d "Density";
    input Temperature T "Temperature";
      input ThermalSeparation.Media.Types.FixedPhase phase=
                             0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output SpecificEnthalpy h "specific enthalpy";
  algorithm
    h :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.h_dT(
        d,
        T,
        phase);
  end specificEnthalpy_dT;

  function specificEnthalpy_ps "Return specific enthalpy from p and s"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
  input ThermalSeparation.Media.Types.FixedPhase phase=
                         0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output SpecificEnthalpy h "specific enthalpy";
  algorithm
    h :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.h_ps(
        p,
        s,
        phase);
  end specificEnthalpy_ps;

  function temperature_ps "Return temperature from p and s"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
      input ThermalSeparation.Media.Types.FixedPhase phase=
                             0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output Temperature T "Temperature";
  algorithm
    T :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.T_ps(
        p,
        s,
        phase);
  end temperature_ps;

  function density_ps "Return density from p and s"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
     input ThermalSeparation.Media.Types.FixedPhase phase=
                            0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output Density d "Density";
  algorithm
   // d := density_psX(p, s, fill(0,0), phase);
     d :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.rho_ps(
        p,
        s,
        phase);
  end density_ps;

  function specificEnthalpy_pT "Return specific enthalpy from p and T"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
      input ThermalSeparation.Media.Types.FixedPhase phase=
                             0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output SpecificEnthalpy h "specific enthalpy";
  algorithm
   h :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.h_pT(
      p, T);
  end specificEnthalpy_pT;

  replaceable function density_pT "Return density from p and T"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
      input ThermalSeparation.Media.Types.FixedPhase phase=
                             0
      "2 for two-phase, 1 for one-phase, 0 if not known";
   output Density d "Density";
  algorithm
    d :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.rho_pT(
      p, T);
  end density_pT;

  function setState_pTX
    "Return thermodynamic state as function of p, T and composition X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:]=reference_X "Mass fractions";
   input ThermalSeparation.Media.Types.FixedPhase phase=
                          0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state "thermodynamic state record";
  algorithm
    state := ThermodynamicState(
      d=density_pT(p,T),
      T=T,
      phase=1,
      h=specificEnthalpy_pT(p,T),
      p=p);
    annotation(Documentation(info="<html></html>"));
  end setState_pTX;

  function setState_phX
    "Return thermodynamic state as function of p, h and composition X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input MassFraction X[:]=reference_X "Mass fractions";
    input ThermalSeparation.Media.Types.FixedPhase phase=
                           0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state "thermodynamic state record";
  algorithm
    state := ThermodynamicState(
        d=density_ph(p, h),
        T=temperature_ph(p, h),
        phase=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.phase_ph(
        p, h),
        h=h,
        p=p);
    annotation(Documentation(info="<html></html>"));
  end setState_phX;

  function setState_psX
    "Return thermodynamic state as function of p, s and composition X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input MassFraction X[:]=reference_X "Mass fractions";
    input ThermalSeparation.Media.Types.FixedPhase phase=
                           0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state "thermodynamic state record";
  algorithm
    state := ThermodynamicState(
        d=density_ps(p, s),
        T=temperature_ps(p, s),
        phase=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.phase_ps(
        p, s),
        h=specificEnthalpy_ps(p, s),
        p=p);

    annotation(Documentation(info="<html></html>"));
  end setState_psX;

  function setState_dTX
    "Return thermodynamic state as function of d, T and composition X or Xi"
    extends Modelica.Icons.Function;
    input Density d "density";
    input Temperature T "Temperature";
    input MassFraction X[:]=reference_X "Mass fractions";
    input ThermalSeparation.Media.Types.FixedPhase phase=
                           0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state "thermodynamic state record";
  algorithm
    state := ThermodynamicState(
        d=d,
        T=T,
        phase=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.phase_dT(
        d, T),
        h=specificEnthalpy_dT(d, T),
        p=pressure_dT(d, T));
  end setState_dTX;

  function dynamicViscosity "Return dynamic viscosity"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output DynamicViscosity eta "Dynamic viscosity";
  algorithm
    eta :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.dynamicViscosity(
        state.d,
        state.T,
        state.p,
        state.phase);
  end dynamicViscosity;

  function thermalConductivity "Return thermal conductivity"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output ThermalConductivity lambda "Thermal conductivity";
  algorithm
    lambda :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.thermalConductivity(
        state.d,
        state.T,
        state.p,
        state.phase);
  end thermalConductivity;

  replaceable function prandtlNumber "Return the Prandtl number"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output PrandtlNumber Pr "Prandtl number";
  algorithm
    Pr := dynamicViscosity(state)*specificHeatCapacityCp(state)/thermalConductivity(
      state);
  end prandtlNumber;

  function pressure "Return pressure"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output AbsolutePressure p "Pressure";
  algorithm
    p := state.p;
  end pressure;

  function temperature "Return temperature"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output Temperature T "Temperature";
  algorithm
    T := state.T;
  end temperature;

  function density "Return density"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output Density d "Density";
  algorithm
    d := state.d;
  end density;

  function specificEnthalpy "Return specific enthalpy"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificEnthalpy h "Specific enthalpy";
  algorithm
    h := state.h;
  end specificEnthalpy;

  function specificInternalEnergy "Return specific internal energy"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificEnergy u "Specific internal energy";
  algorithm
    u := state.h  - state.p/state.d;
  end specificInternalEnergy;

  function specificEntropy "Return specific entropy"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificEntropy s "Specific entropy";
  algorithm
    if dT_explicit then
      s :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.s_dT(
          state.d,
          state.T,
          state.phase);
    elseif pT_explicit then
      s :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.s_pT(
        state.p, state.T);
    else
      s :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.s_ph(
          state.p,
          state.h,
          state.phase);
    end if;
  end specificEntropy;

  function specificGibbsEnergy "Return specific Gibbs energy"
   extends Modelica.Icons.Function;
   input ThermodynamicState state "thermodynamic state record";
   output SpecificEnergy g "Specific Gibbs energy";
  algorithm
   g := state.h - state.T*specificEntropy(state);
  end specificGibbsEnergy;

  function specificHelmholtzEnergy "Return specific Helmholtz energy"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificEnergy f "Specific Helmholtz energy";
  algorithm
    f := state.h - state.p/state.d - state.T*specificEntropy(state);
  end specificHelmholtzEnergy;

  function specificHeatCapacityCp
    "Return specific heat capacity at constant pressure"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificHeatCapacity cp
      "Specific heat capacity at constant pressure";
  algorithm
    if dT_explicit then
      cp :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.cp_dT(
          state.d,
          state.T,
          state.phase);
    elseif pT_explicit then
      cp :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.cp_pT(
        state.p, state.T);
    else
      cp :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.cp_ph(
          state.p,
          state.h,
          state.phase);
    end if;
      annotation (Documentation(info="<html>
                                <p>In the two phase region this function returns the interpolated heat capacity between the
                                liquid and vapour state heat capacities.</p>
                          <html>"));
  end specificHeatCapacityCp;

  function heatCapacity_cp = specificHeatCapacityCp "alias for deprecated name";

  function specificHeatCapacityCv
    "Return specific heat capacity at constant volume"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificHeatCapacity cv "Specific heat capacity at constant volume";
  algorithm
    if dT_explicit then
      cv :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.cv_dT(
          state.d,
          state.T,
          state.phase);
    elseif pT_explicit then
      cv :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.cv_pT(
        state.p, state.T);
    else
      cv :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.cv_ph(
          state.p,
          state.h,
          state.phase);
    end if;
  end specificHeatCapacityCv;

  function heatCapacity_cv = specificHeatCapacityCv "alias for deprecated name";

  function isentropicExponent "Return isentropic exponent"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output IsentropicExponent gamma "Isentropic exponent";
  algorithm
    if dT_explicit then
      gamma :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.isentropicExponent_dT(
          state.d,
          state.T,
          state.phase);
    elseif pT_explicit then
      gamma :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.isentropicExponent_pT(
        state.p, state.T);
    else
      gamma :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.isentropicExponent_ph(
          state.p,
          state.h,
          state.phase);
    end if;
  end isentropicExponent;

  function isentropicEnthalpy "Return isentropic enthalpy"
    extends Modelica.Icons.Function;
    input AbsolutePressure p_downstream "downstream pressure";
    input ThermodynamicState refState "reference state for entropy";
    output SpecificEnthalpy h_is "Isentropic enthalpy";
  algorithm
    h_is :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.isentropicEnthalpy(
        p_downstream,
        specificEntropy(refState),
        0);
  end isentropicEnthalpy;

  function velocityOfSound "Return velocity of sound"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output VelocityOfSound a "Velocity of sound";
  algorithm
    if dT_explicit then
      a :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.velocityOfSound_dT(
          state.d,
          state.T,
          state.phase);
    elseif pT_explicit then
      a :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.velocityOfSound_pT(
        state.p, state.T);
    else
      a :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.velocityOfSound_ph(
          state.p,
          state.h,
          state.phase);
    end if;
  end velocityOfSound;

  function isobaricExpansionCoefficient
    "Return overall the isobaric expansion coefficient beta"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output IsobaricExpansionCoefficient beta "Isobaric expansion coefficient";
  algorithm
    if dT_explicit then
      beta :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.beta_dT(
          state.d,
          state.T,
          state.phase);
    elseif pT_explicit then
      beta :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.beta_pT(
        state.p, state.T);
    else
      beta :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.beta_ph(
          state.p,
          state.h,
          state.phase);
    end if;
  end isobaricExpansionCoefficient;

  function beta = isobaricExpansionCoefficient
    "alias for isobaricExpansionCoefficient for user convenience";

  function isothermalCompressibility
    "Return overall the isothermal compressibility factor"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SI.IsothermalCompressibility kappa "Isothermal compressibility";
  algorithm
    if dT_explicit then
      kappa :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.kappa_dT(
          state.d,
          state.T,
          state.phase);
    elseif pT_explicit then
      kappa :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.kappa_pT(
        state.p, state.T);
    else
      kappa :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.kappa_ph(
          state.p,
          state.h,
          state.phase);
    end if;
  end isothermalCompressibility;

  function kappa = isothermalCompressibility
    "alias of isothermalCompressibility for user convenience";

  // explicit derivative functions for finite element models
  function density_derp_h
    "Return density derivative wrt pressure at const specific enthalpy"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output DerDensityByPressure ddph "Density derivative wrt pressure";
  algorithm
    ddph :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.ddph(
        state.p,
        state.h,
        state.phase);
    annotation(Documentation(info="<html></html>"));
  end density_derp_h;

  function density_derh_p
    "Return density derivative wrt specific enthalpy at constant pressure"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output DerDensityByEnthalpy ddhp "Density derivative wrt specific enthalpy";
  algorithm
    ddhp :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.ddhp(
        state.p,
        state.h,
        state.phase);
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

  function molarMass "Return the molar mass of the medium"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output MolarMass MM "Mixture molar mass";
  algorithm
    MM := fluidConstants[1].molarMass;
  end molarMass;

  replaceable function specificEnthalpy_pTX
    "Return specific enthalpy from p, T, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[nX] "Mass fractions";
    input ThermalSeparation.Media.Types.FixedPhase phase=
                           0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output SpecificEnthalpy h "Specific enthalpy";
  algorithm
    h := specificEnthalpy(setState_pTX(p,T,X,phase));
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
    input MassFraction X[nX] "Mass fractions";
    input ThermalSeparation.Media.Types.FixedPhase phase=
                           0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output Temperature T "Temperature";
  algorithm
    T := temperature(setState_phX(p,h,X,phase));
    annotation(Documentation(info="<html></html>"));
  end temperature_phX;

  replaceable function density_phX "Return density from p, h, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input MassFraction X[nX] "Mass fractions";
      input ThermalSeparation.Media.Types.FixedPhase phase=
                             0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output Density d "Density";
  algorithm
    d := density(setState_phX(p,h,X,phase));
    annotation(Documentation(info="<html></html>"));
  end density_phX;

  replaceable function temperature_psX
    "Return temperature from p,s, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input MassFraction X[nX] "Mass fractions";
    input ThermalSeparation.Media.Types.FixedPhase phase=
                           0 "2 for two-phase, 1 for one-phase, 0 if not known";
   output Temperature T "Temperature";
  algorithm
    T := temperature(setState_psX(p,s,X,phase));
    annotation(Documentation(info="<html></html>"));
  end temperature_psX;

  replaceable function density_psX "Return density from p, s, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input MassFraction X[nX] "Mass fractions";
    input ThermalSeparation.Media.Types.FixedPhase phase=
                           0 "2 for two-phase, 1 for one-phase, 0 if not known";
   output Density d "Density";
  algorithm
    d := density(setState_psX(p,s,X,phase));
    annotation(Documentation(info="<html></html>"));
  end density_psX;

  replaceable function specificEnthalpy_psX
    "Return specific enthalpy from p, s, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input MassFraction X[nX] "Mass fractions";
    input ThermalSeparation.Media.Types.FixedPhase phase=
                           0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output SpecificEnthalpy h "Specific enthalpy";
  algorithm
    h := specificEnthalpy(setState_psX(p,s,X,phase));
    annotation(Documentation(info="<html></html>"));
  end specificEnthalpy_psX;

   function bubbleEnthalpy "Return bubble point specific enthalpy"
      extends Modelica.Icons.Function;
      input SaturationProperties sat "saturation property record";
      output SI.SpecificEnthalpy hl "boiling curve specific enthalpy";
   algorithm
     hl :=
       ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.hl_p(
       sat.psat);
    annotation(Documentation(info="<html></html>"));
   end bubbleEnthalpy;

    function dewEnthalpy "Return dew point specific enthalpy"
      extends Modelica.Icons.Function;
      input SaturationProperties sat "saturation property record";
      output SI.SpecificEnthalpy hv "dew curve specific enthalpy";

    algorithm
      hv :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.hv_p(
        sat.psat);

    end dewEnthalpy;

    function bubbleEntropy "Return bubble point specific entropy"
    extends Modelica.Icons.Function;
    input SaturationProperties sat "saturation property record";
    output SI.SpecificEntropy sl "boiling curve specific entropy";
    algorithm
      sl :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.sl_p(
        sat.psat);

    annotation(Documentation(info="<html></html>"));
    end bubbleEntropy;

    replaceable partial function dewEntropy "Return dew point specific entropy"
    extends Modelica.Icons.Function;
    input SaturationProperties sat "saturation property record";
    output SI.SpecificEntropy sv "dew curve specific entropy";
    algorithm
      sv :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.sv_p(
        sat.psat);
    annotation(Documentation(info="<html></html>"));
    end dewEntropy;

    function bubbleDensity "Return bubble point density"
      extends Modelica.Icons.Function;
      input SaturationProperties sat "saturation property record";
      output Density dl "boiling curve density";
    algorithm
      if ph_explicit or pT_explicit then
        dl :=
          ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.rhol_p(
          sat.psat);
      else
        dl :=
          ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.rhol_T(
          sat.Tsat);
      end if;
    annotation(Documentation(info="<html></html>"));
    end bubbleDensity;

    function dewDensity "Return dew point density"
      extends Modelica.Icons.Function;
      input SaturationProperties sat "saturation property record";
      output Density dv "dew curve density";

    algorithm
      if ph_explicit or pT_explicit then
        dv :=
          ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.rhov_p(
          sat.psat);
      else
        dv :=
          ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.rhov_T(
          sat.Tsat);
      end if;

    end dewDensity;

 function saturationPressure "Return saturation pressure"
      extends Modelica.Icons.Function;
      input Temperature T "temperature";
      output AbsolutePressure p "saturation pressure";
 algorithm
   p :=
     ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Basic.psat(
     T);
    annotation(Documentation(info="<html></html>"));
 end saturationPressure;

    function saturationTemperature "Return saturation temperature"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "pressure";
      output Temperature T "saturation temperature";
    algorithm
      T :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Basic.tsat(
        p);
    annotation(Documentation(info="<html></html>"));
    end saturationTemperature;

    replaceable function saturationPressure_sat "Return saturation temperature"
      extends Modelica.Icons.Function;
      input SaturationProperties sat "saturation property record";
      output AbsolutePressure p "saturation pressure";
    algorithm
      p := sat.psat;
    annotation(Documentation(info="<html></html>"));
    end saturationPressure_sat;

    replaceable function saturationTemperature_sat
    "Return saturation temperature"
      extends Modelica.Icons.Function;
      input SaturationProperties sat "saturation property record";
      output Temperature T "saturation temperature";
    algorithm
      T := sat.Tsat;
    annotation(Documentation(info="<html></html>"));
    end saturationTemperature_sat;

    replaceable partial function saturationTemperature_derp
    "Return derivative of saturation temperature w.r.t. pressure"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "pressure";
      output Real dTp "derivative of saturation temperature w.r.t. pressure";
    algorithm
      dTp :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Basic.dtsatofp(
        p);
    annotation(Documentation(info="<html></html>"));
    end saturationTemperature_derp;

    replaceable function saturationTemperature_derp_sat
    "Return derivative of saturation temperature w.r.t. pressure"
      extends Modelica.Icons.Function;
      input SaturationProperties sat "saturation property record";
      output Real dTp "derivative of saturation temperature w.r.t. pressure";
    algorithm
      dTp := saturationTemperature_derp(sat.psat);
    annotation(Documentation(info="<html></html>"));
    end saturationTemperature_derp_sat;

  replaceable function setState_px
    "Return thermodynamic state from pressure and vapour quality"
    input AbsolutePressure p "Pressure";
    input MassFraction x "Vapour quality";
    output ThermodynamicState state "Thermodynamic state record";
  algorithm
    state := setState_ph(
      p,
      (1-x)*bubbleEnthalpy(setSat_p(p)) + x*dewEnthalpy(setSat_p(p)),
      2);
    annotation(Documentation(info="<html></html>"));
  end setState_px;

  replaceable function setState_Tx
    "Return thermodynamic state from temperature and vapour quality"
    input Temperature T "Temperature";
    input MassFraction x "Vapour quality";
    output ThermodynamicState state "thermodynamic state record";
  algorithm
    state := setState_ph(
      saturationPressure_sat(setSat_T(T)),
      (1-x)*bubbleEnthalpy(setSat_T(T)) + x*dewEnthalpy(setSat_T(T)),
      2);
    annotation(Documentation(info="<html></html>"));
  end setState_Tx;

  replaceable function vapourQuality "Return vapour quality"
    input ThermodynamicState state "Thermodynamic state record";
    output MassFraction x "Vapour quality";
  protected
    constant SpecificEnthalpy eps = 1e-8;
  algorithm
    x := min(max(
      (specificEnthalpy(state)-bubbleEnthalpy(setSat_p(pressure(state)))) /
      (dewEnthalpy(setSat_p(pressure(state))) - bubbleEnthalpy(setSat_p(pressure(state))) + eps),
      0),1);
    annotation(Documentation(info="<html></html>"));
  end vapourQuality;

       function setDewState "Return the thermodynamic state on the dew line"
    extends Modelica.Icons.Function;
    input SaturationProperties sat "saturation point";
    input ThermalSeparation.Media.Types.FixedPhase phase(
      min=1,
      max=2) = 1 "phase: default is one phase";
    output ThermodynamicState state "complete thermodynamic state info";

       algorithm
    state := ThermodynamicState(
       phase = phase,
       p = sat.psat,
       T = sat.Tsat,
       h = dewEnthalpy(sat),
       d = dewDensity(sat));

       end setDewState;

  function setBubbleState "Return the thermodynamic state on the bubble line"
    extends Modelica.Icons.Function;
    input SaturationProperties sat "saturation point";
    input ThermalSeparation.Media.Types.FixedPhase phase(
                           min = 1, max = 2) =  1 "phase: default is one phase";
    output ThermodynamicState state "complete thermodynamic state info";
  algorithm
    state := ThermodynamicState(
       phase = phase,
       p = sat.psat,
       T = sat.Tsat,
       h = bubbleEnthalpy(sat),
       d = bubbleDensity(sat));
    annotation(Documentation(info="<html></html>"));
  end setBubbleState;

  replaceable function setSat_T
    "Return saturation property record from temperature"
    extends Modelica.Icons.Function;
    input Temperature T "temperature";
    output SaturationProperties sat "saturation property record";
  algorithm
    sat.Tsat := T;
    sat.psat := saturationPressure(T);
    annotation(Documentation(info="<html></html>"));
  end setSat_T;

  replaceable function setSat_p
    "Return saturation property record from pressure"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "pressure";
    output SaturationProperties sat "saturation property record";
  algorithm
    sat.psat := p;
    sat.Tsat := saturationTemperature(p);
    annotation(Documentation(info="<html></html>"));
  end setSat_p;

  function surfaceTension
    "Return surface tension sigma in the two phase region"
    extends Modelica.Icons.Function;
    input SaturationProperties sat "saturation property record";
    output SurfaceTension sigma "Surface tension sigma in the two phase region";
  algorithm
    sigma :=
      ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.surfaceTension(
      sat.Tsat);
  end surfaceTension;

    function dBubbleDensity_dPressure "Return bubble point density derivative"
      extends Modelica.Icons.Function;
      input SaturationProperties sat "saturation property record";
      output DerDensityByPressure ddldp "boiling curve density derivative";
    algorithm
      ddldp :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.drhol_dp(
        sat.psat);
    annotation(Documentation(info="<html></html>"));
    end dBubbleDensity_dPressure;

    function dDewDensity_dPressure "Return dew point density derivative"
      extends Modelica.Icons.Function;
      input SaturationProperties sat "saturation property record";
      output DerDensityByPressure ddvdp "saturated steam density derivative";
    algorithm
      ddvdp :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.drhov_dp(
        sat.psat);
    annotation(Documentation(info="<html></html>"));
    end dDewDensity_dPressure;

    function dBubbleEnthalpy_dPressure
    "Return bubble point specific enthalpy derivative"
      extends Modelica.Icons.Function;
      input SaturationProperties sat "saturation property record";
      output DerEnthalpyByPressure dhldp
      "boiling curve specific enthalpy derivative";

    algorithm
      dhldp :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.dhl_dp(
        sat.psat);
    annotation(Documentation(info="<html></html>"));
    end dBubbleEnthalpy_dPressure;

    function dDewEnthalpy_dPressure
    "Return dew point specific enthalpy derivative"
      extends Modelica.Icons.Function;

      input SaturationProperties sat "saturation property record";
      output DerEnthalpyByPressure dhvdp
      "saturated steam specific enthalpy derivative";

    algorithm
      dhvdp :=
        ThermalSeparation.Media.WaterBasedLiquid.BaseClasses.IF97_Utilities.BaseIF97.Regions.dhv_dp(
        sat.psat);
    annotation(Documentation(info="<html></html>"));
    end dDewEnthalpy_dPressure;

  function rho_boiling
    input Integer i;
    output SI.Density rho_boiling "density at boiling temperatur in kg/m3";
   //constant values at 1 bar; O2, N2, SO2, H2O: Refprop, CO2: Refprop extrapolated, HF, HCl: Linde
  protected
    SI.Density rho_boiling_const[:] = {807, 1142, 1300, 959};
  algorithm
    rho_boiling:=rho_boiling_const[i];
  end rho_boiling;

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

      extends Modelica.Icons.Library;
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

      extends Modelica.Icons.Library;
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

      extends Modelica.Icons.Library;
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

      extends Modelica.Icons.Library;
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

      extends Modelica.Icons.Library;
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
end PartialMediumReaction;
