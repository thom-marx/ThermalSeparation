within ThermalSeparation.Media.IdealGasMixtures.BaseClasses;
partial package BaseIdealGasMixture
    extends ThermalSeparation.Media.BaseMediumVapour(MMX=data[:].MM, delta_hv_medium=false);
  import      Modelica.Units.SI;
  extends Modelica.Icons.Package;
  import Modelica.Math;

 constant
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord[
                                                                        nX] data;

  // Constants to be set in Medium
  constant String mediumName = "unusablePartialMedium" "Name of the medium";
  constant String extraPropertiesNames[:]=fill("", 0)
    "Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused";
  constant Boolean singleState= false
    "= true, if u and d are not a function of pressure";
  constant Boolean reducedX=false
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
  constant SpecificEnthalpy h_default = ThermalSeparation.Media.IdealGasMixtures.BaseClasses.BaseIdealGasMixture.specificEnthalpy_pTX(
                                                             p_default, T_default, X_default)
    "Default value for specific enthalpy of medium (for initialization)";
  constant MassFraction X_default[nX]=reference_X
    "Default value for mass fractions of medium (for initialization)";

  constant Integer nS(min=1) "Number of substances"  annotation(Evaluate=true);
  constant Integer nX= nS
    "Number of mass fractions (= 0, if only one substance)" annotation(Evaluate=true);
  constant Integer nXi=if fixedX then 0 else if reducedX then nS - 1 else nX
    "Number of structurally independent mass fractions (see docu for details)"
   annotation(Evaluate=true);

  final constant Integer nC=size(extraPropertiesNames, 1)
    "Number of extra (outside of standard mass-balance) transported properties"
   annotation(Evaluate=true);

   constant SI.Temperature Tcrit[nX]=critTemp();
   constant SI.Pressure pcrit[nX]=critPressure();
   constant SI.MolarVolume Vcrit[nX]=critVolume();
   constant Real omega[nX]=acentricFac();
   constant DipoleMoment mu[nX]=dipoleMoment();
constant Real R_const[nX];

      constant Boolean excludeEnthalpyOfFormation=true
    "If true, enthalpy of formation Hf is not included in specific enthalpy h";
 constant Choices.ReferenceEnthalpy.Temp referenceChoice=Choices.
     ReferenceEnthalpy.ZeroAt0K "Choice of reference enthalpy";
 constant SpecificEnthalpy h_offset=0.0
    "User defined offset for reference enthalpy, if referenceChoice = UserDefined";

   replaceable record FluidConstants "extended fluid constants"
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
     SI.ElectricDipoleMomentOfMolecule dipoleMoment
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
     Boolean hasAcentricFactor=false "true if Pitzer accentric factor is known";
     SpecificEnthalpy HCRIT0=0.0
      "Critical specific enthalpy of the fundamental equation";
     SpecificEntropy SCRIT0=0.0
      "Critical specific entropy of the fundamental equation";
     SpecificEnthalpy deltah=0.0
      "Difference between specific enthalpy model (h_m) and f.eq. (h_f) (h_m - h_f)";
     SpecificEntropy deltas=0.0
      "Difference between specific enthalpy model (s_m) and f.eq. (s_f) (s_m - s_f)";
   end FluidConstants;
 constant FluidConstants[nS] fluidConstants "constant data for the fluid";

  redeclare replaceable model extends BaseProperties
    "Base properties (p, d, T, h, u, R, MM and, if applicable, X) of a medium"

    MassFraction[nXi] Xi(start=reference_X[1:nXi])
      "Structurally independent mass fractions" annotation (Hide=true);
    SpecificHeatCapacity R "Gas constant (of mixture if applicable)";
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

     /*** calculation of eta ***/
  protected
      Real fi[nX, nX];
      DynamicViscosity[nX] eta_pure;
    parameter Real Const1_SI=40.785*10^(-9.5)
      "Constant in formula for eta converted to SI units";
    parameter Real Const2_SI=131.3/1000.0
      "Constant in formula for mur converted to SI units";
    Real mur[nX] "Dimensionless dipole moment of gas molecule";
    Real Fc[nX]=ones(nX) - 0.2756*omega + 0.059035*mur.^4
      "Factor to account for molecular shape and polarities of gas";
    Real Tstar[nX] "Dimensionless temperature defined by equation below";
    Real Ov[nX] "Viscosity collision integral for the gas";

  equation
    /*** calculation of eta ***/
   for i in 1:nX loop
        mur[i] = Const2_SI*mu[i]/sqrt(Vcrit[i]*Tcrit[i]);
              Tstar[i] =  1.2593*state.T/Tcrit[i];
                Ov[i] =  1.16145*Tstar[i]^(-0.14874) + 0.52487*exp(-0.7732*Tstar[i]) + 2.16178*exp(-2.43787*Tstar[i]);
        eta_pure[i] = Const1_SI*Fc[i]*sqrt(MMX[i]*state.T)/(Vcrit[i]^(2/3)*Ov[i]);
    end for;
      for i in 1:nX loop
      for j in 1:nX loop
        if i==1 then
          fi[i,j] =  (1 + (eta_pure[i]/eta_pure[j])^(1/2)*(MMX[j]/MMX[i])^(1/4))^2/(8*(1 + MMX[i]/MMX[j]))^(1/2);
        elseif j<i then
            fi[i,j] =  eta_pure[i]/eta_pure[j]*MMX[j]/MMX[i]*fi[j,i];
          else
            fi[i,j] =  (1 + (eta_pure[i]/eta_pure[j])^(1/2)*(MMX[j]/MMX[i])^(1/4))^2/(8*(1 + MMX[i]/MMX[j]))^(1/2);
        end if;
      end for;
    end for;
      eta =sum(x[i]*eta_pure[i]/max(1e-5, sum(x[j]*fi[i, j] for j in 1:nX)) for i in 1:nX);

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
           assert(X[i] >= -1.e-5 and X[i] <= 1 + 1.e-5, "Mass fraction X is not in the range 0..1");
         end for;
       end if;
     end if;

     assert(p >= 0.0, "Pressure (= " + String(p) + " Pa) of medium \"" +
       mediumName + "\" is negative\n(Temperature = " + String(T) + " K)");

           lambda = thermalConductivity(state);
      cp=specificHeatCapacityCp(state);

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

 function gasConstant
    "Return the gas constant of the mixture (also for liquids)"
     extends Modelica.Icons.Function;
     input ThermodynamicState state "thermodynamic state";
     output SI.SpecificHeatCapacity R "mixture gas constant";
 algorithm
     R := R_const*state.X;
 end gasConstant;

   function moleToMassFractions "Return mass fractions X from mole fractions"
     extends Modelica.Icons.Function;
     input SI.MoleFraction moleFractions[:] "Mole fractions of mixture";
     input MolarMass[:] MMX "molar masses of components";
     output SI.MassFraction X[size(moleFractions, 1)]
      "Mass fractions of gas mixture";
  protected
     MolarMass Mmix =  moleFractions*MMX "molar mass of mixture";
   algorithm
     for i in 1:size(moleFractions, 1) loop
       X[i] := moleFractions[i]*MMX[i] /Mmix;
     end for;
   end moleToMassFractions;

   function massToMoleFractions "Return mole fractions from mass fractions X"
     extends Modelica.Icons.Function;
     input SI.MassFraction X[:] "Mass fractions of mixture";
     input SI.MolarMass[:] MMX "molar masses of components";
     output SI.MoleFraction moleFractions[size(X, 1)]
      "Mole fractions of gas mixture";
  protected
     Real invMMX[size(X, 1)] "inverses of molar weights";
     SI.MolarMass Mmix "molar mass of mixture";
   algorithm
     for i in 1:size(X, 1) loop
       invMMX[i] := 1/MMX[i];
     end for;
     Mmix := 1/(X*invMMX);
     for i in 1:size(X, 1) loop
       moleFractions[i] := Mmix*X[i]/MMX[i];
     end for;
   end massToMoleFractions;

         function setState_pTX "Return thermodynamic state as function of p, T and composition X"
           extends Modelica.Icons.Function;
           input AbsolutePressure p "Pressure";
           input Temperature T "Temperature";
           input MassFraction X[:]=reference_X "Mass fractions";
           output ThermodynamicState state;
         algorithm
           state := if size(X, 1) == nX then ThermodynamicState(
               p=p,
               T=T,
               X=X) else
            ThermodynamicState(p=p,T=T, X=cat(1,X,{1-sum(X)}));
         end setState_pTX;

       function setState_phX "Return thermodynamic state as function of p, h and composition X"
         extends Modelica.Icons.Function;
         input AbsolutePressure p "Pressure";
         input SpecificEnthalpy h "Specific enthalpy";
         input MassFraction X[:]=reference_X "Mass fractions";
         output ThermodynamicState state;
       algorithm
         state := if size(X, 1) == nX then ThermodynamicState(
             p=p,
             T=T_hX(h, X),
             X=X) else
            ThermodynamicState(p=p,T=T_hX(h,X), X=cat(1,X,{1-sum(X)}));
       end setState_phX;

   function setState_psX
    "Return thermodynamic state as function of p, s and composition X"
     extends Modelica.Icons.Function;
     input AbsolutePressure p "Pressure";
     input SpecificEntropy s "Specific entropy";
     input MassFraction X[:]=reference_X "Mass fractions";
     output ThermodynamicState state;
   algorithm
     state := if size(X,1) == nX then ThermodynamicState(p=p,T=T_psX(p,s,X),X=X) else
            ThermodynamicState(p=p,T=T_psX(p,s,X), X=cat(1,X,{1-sum(X)}));
   end setState_psX;

   function setState_dTX
    "Return thermodynamic state as function of d, T and composition X"
     extends Modelica.Icons.Function;
     input Density d "density";
     input Temperature T "Temperature";
     input MassFraction X[nX]=reference_X "Mass fractions";
     output ThermodynamicState state;

   algorithm
        state := if size(X,1) == nX then ThermodynamicState(p=d*(R_const*X)*T,T=T,X=X) else
            ThermodynamicState(p=d*(R_const*cat(1,X,{1-sum(X)}))*T,T=T, X=cat(1,X,{1-sum(X)}));

   end setState_dTX;

   function pressure "Return pressure of ideal gas"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output AbsolutePressure p "Pressure";
   algorithm
    p := state.p;
   end pressure;

   function temperature "Return temperature of ideal gas"
      extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output Temperature T "Temperature";
   algorithm
    T := state.T;
   end temperature;

   function density "Return density of ideal gas"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output Density d "Density";
   algorithm
    d := state.p/((state.X*R_const)*state.T);
   end density;

 function specificEnthalpy "Return specific enthalpy"
    extends Modelica.Icons.Function;
   input ThermodynamicState state "thermodynamic state record";
   output SpecificEnthalpy h "Specific enthalpy";
 algorithm
   h := h_TX(state.T,state.X);
 end specificEnthalpy;

 function specificInternalEnergy "Return specific internal energy"
    extends Modelica.Icons.Function;
  input ThermodynamicState state "thermodynamic state record";
  output SpecificEnergy u "Specific internal energy";
 algorithm
  u := h_TX(state.T,state.X) - gasConstant(state)*state.T;
 end specificInternalEnergy;

 function specificEntropy "Return specific entropy"
   extends Modelica.Icons.Function;
  input ThermodynamicState state "thermodynamic state record";
  output SpecificEntropy s "Specific entropy";
  protected
  Real[nX] Y(unit="mol/mol")=massToMoleFractions(state.X, MMX)
      "Molar fractions";
 algorithm
  s := s_TX(state.T, state.X) - sum(state.X[i]*Modelica.Constants.R/MMX[i]*
    Modelica.Math.log(Y[i]*state.p/
    reference_p) for i in 1:nX);
    //  s := s_TX(state.T, state.X) - gasConstant(state)*Modelica.Math.log(state.p/reference_p)
    //    + MixEntropy(massToMoleFractions(state.X,data.MM));
 end specificEntropy;

function specificGibbsEnergy "Return specific Gibbs energy"
    extends Modelica.Icons.Function;
 input ThermodynamicState state "thermodynamic state record";
 output SpecificEnergy g "Specific Gibbs energy";
algorithm
  g := h_TX(state.T,state.X) - state.T*specificEntropy(state);
end specificGibbsEnergy;

 function specificHelmholtzEnergy "Return specific Helmholtz energy"
     extends Modelica.Icons.Function;
   input ThermodynamicState state "thermodynamic state record";
   output SpecificEnergy f "Specific Helmholtz energy";
 algorithm
   f := h_TX(state.T,state.X) - gasConstant(state)*state.T - state.T*specificEntropy(state);
 end specificHelmholtzEnergy;

 function h_TX "Return specific enthalpy"
    import
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices;
    extends Modelica.Icons.Function;
    input SI.Temperature T "Temperature";
    input MassFraction X[:]=reference_X
      "Independent Mass fractions of gas mixture";
    input Boolean exclEnthForm=excludeEnthalpyOfFormation
      "If true, enthalpy of formation Hf is not included in specific enthalpy h";
    input Choices.ReferenceEnthalpy.Temp refChoice=referenceChoice
      "Choice of reference enthalpy";
    input SI.SpecificEnthalpy h_off=h_offset
      "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
    output SI.SpecificEnthalpy h "Specific enthalpy at temperature T";
 algorithm
   h := (if fixedX then reference_X else X)*{
     ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasNasa.h_T(
     data[i],
     T,
     exclEnthForm,
     refChoice,
     h_off) for i in 1:nX};
   annotation(Inline=false,smoothOrder=1);
 end h_TX;

 function h_TX_der "Return specific enthalpy derivative"
    import
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices;
    extends Modelica.Icons.Function;
    input SI.Temperature T "Temperature";
    input MassFraction X[nX] "Independent Mass fractions of gas mixture";
    input Boolean exclEnthForm=excludeEnthalpyOfFormation
      "If true, enthalpy of formation Hf is not included in specific enthalpy h";
    input Choices.ReferenceEnthalpy.Temp refChoice=referenceChoice
      "Choice of reference enthalpy";
    input SI.SpecificEnthalpy h_off=h_offset
      "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
   input Real dT "Temperature derivative";
   input Real dX[nX] "independent mass fraction derivative";
   output Real h_der "Specific enthalpy at temperature T";
 algorithm
   h_der := if fixedX then dT*sum((
     ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.Functions.cp_T(
      data[i], T)*reference_X[i]) for i in 1:nX) else dT*sum((
     ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.Functions.cp_T(
      data[i], T)*X[i]) for i in 1:nX) + sum((
     ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasNasa.h_T(
      data[i], T)*dX[i]) for i in 1:nX);
   annotation (InlineNoEvent=false, Inline = false);
 end h_TX_der;
 // redeclare function extends gasConstant "Return gasConstant"
 // algorithm
 //   R := data.R_s*state.X;
 // end gasConstant;

 function specificHeatCapacityCp
    "Return specific heat capacity at constant pressure"
     extends Modelica.Icons.Function;
   input ThermodynamicState state "thermodynamic state record";
   output SpecificHeatCapacity cp "Specific heat capacity at constant pressure";
 algorithm
    cp := {
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.Functions.cp_T(
      data[i], state.T) for i in 1:nX}*state.X;
 end specificHeatCapacityCp;

 function heatCapacity_cp = specificHeatCapacityCp "alias for deprecated name";
 function specificHeatCapacityCv
    "Return specific heat capacity at constant volume from temperature and gas data"
     extends Modelica.Icons.Function;
   input ThermodynamicState state "thermodynamic state record";
   output SpecificHeatCapacity cv "Specific heat capacity at constant volume";
 algorithm
    cv := {
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.Functions.cp_T(
      data[i], state.T) for i in 1:nX}*state.X - data.R_s*state.X;
 end specificHeatCapacityCv;

 function heatCapacity_cv = specificHeatCapacityCv "alias for deprecated name";
 function MixEntropy "Return mixing entropy of ideal gases / R"
   extends Modelica.Icons.Function;
   input SI.MoleFraction x[:] "mole fraction of mixture";
   output Real smix "mixing entropy contribution, divided by gas constant";
 algorithm
   smix := sum(if x[i] > Modelica.Constants.eps then -x[i]*Modelica.Math.log(x[i]) else
                    x[i] for i in 1:size(x,1));
 end MixEntropy;

 function s_TX
    "Return temperature dependent part of the entropy, expects full entropy vector"
   input Temperature T "temperature";
   input MassFraction[nX] X "mass fraction";
   output SpecificEntropy s "specific entropy";
 algorithm
    s := sum(
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasNasa.s0_T(
      data[i], T)*X[i] for i in 1:size(X, 1));
 end s_TX;

 function isentropicExponent "Return isentropic exponent"
   extends Modelica.Icons.Function;
   input ThermodynamicState state "thermodynamic state record";
   output IsentropicExponent gamma "Isentropic exponent";
 algorithm
   gamma := specificHeatCapacityCp(state)/specificHeatCapacityCv(state);
 end isentropicExponent;

 function velocityOfSound "Return velocity of sound"
   extends Modelica.Icons.Function;
   input ThermodynamicState state "properties at upstream location";
   output VelocityOfSound a "Velocity of sound";
 algorithm
   a := sqrt(gasConstant(state)*state.T*specificHeatCapacityCp(state)/specificHeatCapacityCv(state));
 end velocityOfSound;

 function isentropicEnthalpyApproximation
    "Approximate method of calculating h_is from upstream properties and downstream pressure"
   extends Modelica.Icons.Function;
   input AbsolutePressure p2 "downstream pressure";
   input ThermodynamicState state "thermodynamic state at upstream location";
   output SpecificEnthalpy h_is "isentropic enthalpy";
  protected
   SpecificEnthalpy h "specific enthalpy at upstream location";
   SpecificEnthalpy h_component[nX] "specific enthalpy at upstream location";
   IsentropicExponent gamma =  isentropicExponent(state) "Isentropic exponent";
  protected
   MassFraction[nX] X "complete X-vector";
 algorithm
   X := if reducedX then cat(1,state.X,{1-sum(state.X)}) else state.X;
    h_component := {
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.SingleGasNasa.h_T(
       data[i],
       state.T,
       excludeEnthalpyOfFormation,
       referenceChoice,
       h_offset) for i in 1:nX};
   h :=h_component*X;
   h_is := h + gamma/(gamma - 1.0)*(state.T*gasConstant(state))*
     ((p2/state.p)^((gamma - 1)/gamma) - 1.0);
 end isentropicEnthalpyApproximation;

 function isentropicEnthalpy "Return isentropic enthalpy"
   extends Modelica.Icons.Function;
   input AbsolutePressure p_downstream "downstream pressure";
   input ThermodynamicState refState "reference state for entropy";
   input Boolean exact = false
      "flag wether exact or approximate version should be used";
   output SpecificEnthalpy h_is "Isentropic enthalpy";

 algorithm
   h_is := if exact then specificEnthalpy_psX(p_downstream,specificEntropy(refState),refState.X) else
          isentropicEnthalpyApproximation(p_downstream,refState);
 end isentropicEnthalpy;

  function gasMixtureViscosity
    "Return viscosities of gas mixtures at low pressures (Wilke method)"
    extends Modelica.Icons.Function;
    input MoleFraction[:] yi "Mole fractions";
    input MolarMass[:] M "Mole masses";
    input DynamicViscosity[:] eta "Pure component viscosities";
    output DynamicViscosity etam "Viscosity of the mixture";
  protected
    Real fi[size(yi,1),size(yi,1)];
  algorithm
    for i in 1:size(eta,1) loop
      assert(fluidConstants[i].hasDipoleMoment,"Dipole moment for " + fluidConstants[i].chemicalFormula +
         " not known. Can not compute viscosity.");
      assert(fluidConstants[i].hasCriticalData, "Critical data for "+ fluidConstants[i].chemicalFormula +
         " not known. Can not compute viscosity.");
      for j in 1:size(eta,1) loop
        if i==1 then
          fi[i,j] := (1 + (eta[i]/eta[j])^(1/2)*(M[j]/M[i])^(1/4))^2/(8*(1 + M[i]/M[j]))^(1/2);
        elseif j<i then
            fi[i,j] := eta[i]/eta[j]*M[j]/M[i]*fi[j,i];
          else
            fi[i,j] := (1 + (eta[i]/eta[j])^(1/2)*(M[j]/M[i])^(1/4))^2/(8*(1 + M[i]/M[j]))^(1/2);
        end if;
      end for;
    end for;
    etam := sum(yi[i]*eta[i]/sum(yi[j]*fi[i,j] for j in 1:size(eta,1)) for i in 1:size(eta,1));

   annotation (Documentation(info="<html>
 
<p>
Simplification of the kinetic theory (Chapman and Enskog theory)
approach neglecting the second-order effects.<br>
<br>
This equation has been extensively tested (Amdur and Mason, 1958;
Bromley and Wilke, 1951; Cheung, 1958; Dahler, 1959; Gandhi and Saxena,
1964; Ranz and Brodowsky, 1962; Saxena and Gambhir, 1963a; Strunk, et
al., 1964; Vanderslice, et al. 1962; Wright and Gray, 1962). In most
cases, only nonpolar mixtures were compared, and very good results
obtained. For some systems containing hidrogen as one component, less
satisfactory agreement was noted. Wilke's method predicted mixture
viscosities that were larger than experimental for the H2-N2 system,
but for H2-NH3, it underestimated the viscosities. <br>
Gururaja, et al. (1967) found that this method also overpredicted in
the H2-O2 case but was quite accurate for the H2-CO2 system. <br>
Wilke's approximation has proved reliable even for polar-polar gas
mixtures of aliphatic alcohols (Reid and Belenyessy, 1960). The
principal reservation appears to lie in those cases where Mi&gt;&gt;Mj
and etai&gt;&gt;etaj.<br>
</p>
 
</html>
"));
  end gasMixtureViscosity;

   function dynamicViscosity "Return mixture dynamic viscosity"
     extends Modelica.Icons.Function;
     input ThermodynamicState state "thermodynamic state record";
     output DynamicViscosity eta "Dynamic viscosity";
  protected
     DynamicViscosity[nX] etaX "component dynamic viscosities";
   algorithm
     for i in 1:nX loop
      etaX[i] :=
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.Functions.dynamicViscosityLowPressure(
           state.T,
           fluidConstants[i].criticalTemperature,
           fluidConstants[i].molarMass,
           fluidConstants[i].criticalMolarVolume,
           fluidConstants[i].acentricFactor,
           fluidConstants[i].dipoleMoment);
     end for;
     eta := gasMixtureViscosity(massToMoleFractions(state.X,
                            fluidConstants[:].molarMass),
                fluidConstants[:].molarMass,
                etaX);
   end dynamicViscosity;

 function mixtureViscosityChung
    "Return the viscosity of gas mixtures without access to component viscosities (Chung, et. al. rules)"
 extends Modelica.Icons.Function;
    import      Modelica.Units.SI;
   input Temperature T "Temperature";
   input Temperature[:] Tc "Critical temperatures";
   input MolarVolume[:] Vcrit "Critical volumes (m3/mol)";
   input Real[:] w "Acentric factors";
   input Real[:] mu "Dipole moments (debyes)";
   input MolarMass[:] MolecularWeights "Molecular weights (kg/mol)";
   input MoleFraction[:] y "Molar Fractions";
   input Real[:] kappa =  zeros(nX) "Association Factors";
   output DynamicViscosity etaMixture "Mixture viscosity (Pa.s)";
  protected
 constant Real[size(y,1)] Vc =  Vcrit*1000000 "Critical volumes (cm3/mol)";
 constant Real[size(y,1)] M =  MolecularWeights*1000
      "Molecular weights (g/mol)";
 Integer n = size(y,1) "Number of mixed elements";
 Real sigmam3 "Mixture sigma^3 in Ångström";
 Real sigma[size(y,1),size(y,1)];
 Real edivkm;
 Real edivk[size(y,1),size(y,1)];
 Real Mm;
 Real Mij[size(y,1),size(y,1)];
 Real wm "accentric factor";
 Real wij[size(y,1),size(y,1)];
 Real kappam
      "Correlation for highly polar substances such as alcohols and acids";
 Real kappaij[size(y,1),size(y,1)];
 Real mum;
 Real Vcm;
 Real Tcm;
 Real murm "Dimensionless dipole moment of the mixture";
 Real Fcm "Factor to correct for shape and polarity";
 Real omegav;
 Real Tmstar;
 Real etam "Mixture viscosity in microP";
 algorithm
 //combining rules
 for i in 1:n loop
   for j in 1:n loop
     Mij[i,j] := 2*M[i]*M[j]/(M[i]+M[j]);
     if i==j then
       sigma[i,j] := 0.809*Vc[i]^(1/3);
       edivk[i,j] := Tc[i]/1.2593;
       wij[i,j] := w[i];
       kappaij[i,j] := kappa[i];
     else
       sigma[i,j] := (0.809*Vc[i]^(1/3)*0.809*Vc[j]^(1/3))^(1/2);
       edivk[i,j] := (Tc[i]/1.2593*Tc[j]/1.2593)^(1/2);
       wij[i,j] := (w[i] + w[j])/2;
       kappaij[i,j] := (kappa[i]*kappa[j])^(1/2);
     end if;
   end for;
 end for;
 //mixing rules
 sigmam3 := (sum(sum(y[i]*y[j]*sigma[i,j]^3 for j in 1:n) for i in 1:n));
 //(epsilon/k)m
 edivkm := (sum(sum(y[i]*y[j]*edivk[i,j]*sigma[i,j]^3 for j in 1:n) for i in 1:n))/sigmam3;
 Mm := ((sum(sum(y[i]*y[j]*edivk[i,j]*sigma[i,j]^2*Mij[i,j]^(1/2) for j in 1:n) for i in 1:n))/(edivkm*sigmam3^(2/3)))^2;
 wm := (sum(sum(y[i]*y[j]*wij[i,j]*sigma[i,j]^3 for j in 1:n) for i in 1:n))/sigmam3;
 mum := (sigmam3*(sum(sum(y[i]*y[j]*mu[i]^2*mu[j]^2/sigma[i,j]^3 for j in 1:n) for i in 1:n)))^(1/4);
 Vcm := sigmam3/(0.809)^3;
 Tcm := 1.2593*edivkm;
 murm := 131.3*mum/(Vcm*Tcm)^(1/2);
 kappam := (sigmam3*(sum(sum(y[i]*y[j]*kappaij[i,j] for j in 1:n) for i in 1:n)));
 Fcm := 1 - 0.275*wm + 0.059035*murm^4 + kappam;
 Tmstar := T/edivkm;
 omegav := 1.16145*(Tmstar)^(-0.14874) + 0.52487*Math.exp(-0.77320*Tmstar) + 2.16178*Math.exp(-2.43787*Tmstar);
 etam := 26.69*Fcm*(Mm*T)^(1/2)/(sigmam3^(2/3)*omegav);
 etaMixture := etam*1e7;

 annotation (Documentation(info="<html>
 
<p>
Equation to estimate the viscosity of gas mixtures at low pressures.<br>
It is a simplification of an extension of the rigorous kinetic theory
of Chapman and Enskog to determine the viscosity of multicomponent
mixtures, at low pressures and with a factor to correct for molecule
shape and polarity.
</p>
 
<p>
The input argument Kappa is a special correction for highly polar substances such as
alcohols and acids.<br>
Values of kappa for a few such materials:
</p>
 
<table style=\"text-align: left; width: 302px; height: 200px;\" border=\"1\"
cellspacing=\"0\" cellpadding=\"2\">
<tbody>
<tr>
<td style=\"vertical-align: top;\">Compound <br>
</td>
<td style=\"vertical-align: top; text-align: center;\">Kappa<br>
</td>
<td style=\"vertical-align: top;\">Compound<br>
</td>
<td style=\"vertical-align: top;\">Kappa<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">Methanol<br>
</td>
<td style=\"vertical-align: top;\">0.215<br>
</td>
<td style=\"vertical-align: top;\">n-Pentanol<br>
</td>
<td style=\"vertical-align: top;\">0.122<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">Ethanol<br>
</td>
<td style=\"vertical-align: top;\">0.175<br>
</td>
<td style=\"vertical-align: top;\">n-Hexanol<br>
</td>
<td style=\"vertical-align: top;\">0.114<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">n-Propanol<br>
</td>
<td style=\"vertical-align: top;\">0.143<br>
</td>
<td style=\"vertical-align: top;\">n-Heptanol<br>
</td>
<td style=\"vertical-align: top;\">0.109<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">i-Propanol<br>
</td>
<td style=\"vertical-align: top;\">0.143<br>
</td>
<td style=\"vertical-align: top;\">Acetic Acid<br>
</td>
<td style=\"vertical-align: top;\">0.0916<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">n-Butanol<br>
</td>
<td style=\"vertical-align: top;\">0.132<br>
</td>
<td style=\"vertical-align: top;\">Water<br>
</td>
<td style=\"vertical-align: top;\">0.076<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">i-Butanol<br>
</td>
<td style=\"vertical-align: top;\">0.132</td>
<td style=\"vertical-align: top;\"><br>
</td>
<td style=\"vertical-align: top;\"><br>
</td>
</tr>
</tbody>
</table>
<p>
Chung, et al. (1984) suggest that for other alcohols not shown in the
table:<br>
&nbsp;&nbsp;&nbsp;&nbsp; <br>
&nbsp;&nbsp;&nbsp; kappa = 0.0682 + 4.704*[(number of -OH
groups)]/[molecular weight]<br>
<br>
<span style=\"font-weight: normal;\">S.I. units relation for the
debyes:&nbsp;</span><br>
&nbsp;&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp;
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp;
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp;
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; 1 debye = 3.162e-25 (J.m^3)^(1/2)<br>
</p>
<h4>References</h4>
<p>
[1] THE PROPERTIES OF GASES AND LIQUIDS, Fifth Edition,<br>
&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp; Bruce E. Poling, John M.
Prausnitz, John P. O'Connell.<br>
[2] Chung, T.-H., M. Ajlan, L. L. Lee, and K. E. Starling: Ind. Eng.
Chem. Res., 27: 671 (1988).<br>
[3] Chung, T.-H., L. L. Lee, and K. E. Starling; Ing. Eng. Chem.
Fundam., 23: 3 ()1984).<br>
</p>
</html>
"));
 end mixtureViscosityChung;

  function lowPressureThermalConductivity
    "Return thermal conductivites of low-pressure gas mixtures (Mason and Saxena Modification)"
    extends Modelica.Icons.Function;
    input MoleFraction[:] y
      "Mole fraction of the components in the gass mixture";
    input Temperature T "Temperature";
    input Temperature[:] Tc "Critical temperatures";
    input AbsolutePressure[:] Pc "Critical pressures";
    input MolarMass[:] M "Molecular weights";
    input ThermalConductivity[:] lambda
      "Thermal conductivities of the pure gases";
    output ThermalConductivity lambdam
      "Thermal conductivity of the gas mixture";
  protected
    MolarMass[size(y,1)] gamma;
    Real[size(y,1)] Tr "Reduced temperature";
    Real[size(y,1),size(y,1)] A "Mason and Saxena Modification";
    constant Real epsilon =  1.0 "Numerical constant near unity";
  algorithm
    for i in 1:size(y,1) loop
      gamma[i] := 210*(Tc[i]*M[i]^3/Pc[i]^4)^(1/6);
      Tr[i] := T/Tc[i];
    end for;
    for i in 1:size(y,1) loop
      for j in 1:size(y,1) loop
        A[i,j] := epsilon*(1 + (gamma[j]*(Math.exp(0.0464*Tr[i]) - Math.exp(-0.2412*Tr[i]))/
        (gamma[i]*(Math.exp(0.0464*Tr[j]) - Math.exp(-0.2412*Tr[j]))))^(1/2)*(M[i]/M[j])^(1/4))^2/
        (8*(1 + M[i]/M[j]))^(1/2);
      end for;
    end for;
    lambdam := sum(y[i]*lambda[i]/(sum(y[j]*A[i,j] for j in 1:size(y,1))) for i in 1:size(y,1));

    annotation (Documentation(info="<html>
 
<p>
This function applies the Masson and Saxena modification of the
Wassiljewa Equation for the thermal conductivity for gas mixtures of
n elements at low pressure.
</p>
 
<p>
For nonpolar gas mixtures errors will generally be less than 3 to 4%.
For mixtures of nonpolar-polar and polar-polar gases, errors greater
than 5 to 8% may be expected. For mixtures in which the sizes and
polarities of the constituent molecules are not greatly different, the
thermal conductivity can be estimated satisfactorily by a mole fraction
average of the pure component conductivities.
</p>
 
</html>
"));
  end lowPressureThermalConductivity;

   function thermalConductivity
    "Return thermal conductivity for low pressure gas mixtures"
      extends Modelica.Icons.Function;
     input ThermodynamicState state "thermodynamic state record";
       input Integer method=1
      "method to compute single component thermal conductivity";
     output ThermalConductivity lambda "Thermal conductivity";
  protected
     ThermalConductivity[nX] lambdaX "component thermal conductivities";
     DynamicViscosity[nX] eta "component thermal dynamic viscosities";
     SpecificHeatCapacity[nX] cp "component heat capacity";
   algorithm
     for i in 1:nX loop
   assert(fluidConstants[i].hasCriticalData, "Critical data for "+ fluidConstants[i].chemicalFormula +
      " not known. Can not compute thermal conductivity.");
      eta[i] :=
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.Functions.dynamicViscosityLowPressure(
           state.T,
           fluidConstants[i].criticalTemperature,
           fluidConstants[i].molarMass,
           fluidConstants[i].criticalMolarVolume,
           fluidConstants[i].acentricFactor,
           fluidConstants[i].dipoleMoment);
      cp[i] :=
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.Functions.cp_T(
        data[i], state.T);
      lambdaX[i] :=
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.Functions.thermalConductivityEstimate(
           Cp=cp[i],
           eta=eta[i],
           method=method,
           MM=fluidConstants[i].molarMass);
     end for;
     lambda := lowPressureThermalConductivity(massToMoleFractions(state.X,
                                  fluidConstants[:].molarMass),
                          state.T,
                          fluidConstants[:].criticalTemperature,
                          fluidConstants[:].criticalPressure,
                          fluidConstants[:].molarMass,
                          lambdaX);
   end thermalConductivity;

 function isobaricExpansionCoefficient
    "Return isobaric expansion coefficient beta"
     extends Modelica.Icons.Function;
   input ThermodynamicState state "thermodynamic state record";
   output IsobaricExpansionCoefficient beta "Isobaric expansion coefficient";
 algorithm
   beta := 1/state.T;
 end isobaricExpansionCoefficient;

 function beta = isobaricExpansionCoefficient
    "alias for isobaricExpansionCoefficient for user convenience";
 function isothermalCompressibility "Return isothermal compressibility factor"
   extends Modelica.Icons.Function;
   input ThermodynamicState state "thermodynamic state record";
   output SI.IsothermalCompressibility kappa "Isothermal compressibility";
 algorithm
   kappa := 1.0/state.p;
 end isothermalCompressibility;

 function kappa = isothermalCompressibility
    "alias of isothermalCompressibility for user convenience";
 function density_derp_T
    "Return density derivative by temperature at constant pressure"
     extends Modelica.Icons.Function;
   input ThermodynamicState state "thermodynamic state record";
   output DerDensityByPressure ddpT "Density derivative wrt pressure";
 algorithm
   ddpT := 1/(state.T*gasConstant(state));
 end density_derp_T;

 function density_derT_p
    "Return density derivative by temperature at constant pressure"
     extends Modelica.Icons.Function;
   input ThermodynamicState state "thermodynamic state record";
   output DerDensityByTemperature ddTp "Density derivative wrt temperature";
 algorithm
   ddTp := -state.p/(state.T*state.T*gasConstant(state));
 end density_derT_p;

      function density_derX "Return density derivative by mass fraction"
         extends Modelica.Icons.Function;
       input ThermodynamicState state "thermodynamic state record";
       output Density[nX] dddX "Derivative of density wrt mass fraction";
      algorithm
       dddX := {-state.p/(state.T*gasConstant(state))*molarMass(state)/data[i].MM for
        i in 1:nX};
      end density_derX;

 function molarMass "Return molar mass of mixture"
  extends Modelica.Icons.Function;
  input ThermodynamicState state "thermodynamic state record";
  output MolarMass MM "Mixture molar mass";
 algorithm
  MM := 1/sum(state.X[j]/data[j].MM for j in 1:size(state.X, 1));
 end molarMass;

 function T_hX "Return temperature from specific enthalpy and mass fraction"
   input SpecificEnthalpy h "specific enthalpy";
   input MassFraction[:] X "mass fractions of composition";
   output Temperature T "temperature";
  protected
   MassFraction[nX] Xfull = if size(X,1) == nX then X else cat(1,X,{1-sum(X)});
 package Internal
      "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
   extends ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.OneNonLinearEquation;
   redeclare record extends f_nonlinear_Data
        "Data to be passed to non-linear function"
     extends ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord;
   end f_nonlinear_Data;

   redeclare function extends f_nonlinear
   algorithm
       y := h_TX(x,X);
   end f_nonlinear;

   // Dummy definition has to be added for current Dymola
   redeclare function extends solve
   end solve;
 end Internal;
 algorithm
   T := Internal.solve(h, 200, 6000, 1.0e5, Xfull, data[1]);
 end T_hX;

  function T_psX
    "Return temperature from pressure, specific entropy and mass fraction"
      input AbsolutePressure p "pressure";
      input SpecificEntropy s "specific entropy";
      input MassFraction[:] X "mass fractions of composition";
      output Temperature T "temperature";
  protected
      MassFraction[nX] Xfull = if size(X,1) == nX then X else cat(1,X,{1-sum(X)});
    package Internal
      "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
      extends ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.OneNonLinearEquation;
      redeclare record extends f_nonlinear_Data
        "Data to be passed to non-linear function"
        extends ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord;
      end f_nonlinear_Data;

      redeclare function extends f_nonlinear
        "note that this function always sees the complete mass fraction vector"
      algorithm
        y := s_TX(x,X)- data[:].R*X*Modelica.Math.log(p/reference_p)
          + MixEntropy(massToMoleFractions(X,data[:].MM));
      end f_nonlinear;

      // Dummy definition has to be added for current Dymola
      redeclare function extends solve
      end solve;
    end Internal;
  algorithm
      T := Internal.solve(s, 200, 6000, p, Xfull, data[1]);
  end T_psX;
  //   redeclare function extends specificEnthalpy_psX
  //   protected
  //     Temperature T "temperature";
  //   algorithm
  //     T := temperature_psX(p,s,X);
  //     h := specificEnthalpy_pTX(p,T,X);
  //   end extends;

  //   redeclare function extends density_phX
  //     "Compute density from pressure, specific enthalpy and mass fraction"
  //     protected
  //     Temperature T "temperature";
  //     SpecificHeatCapacity R "gas constant";
  //   algorithm
  //     T := temperature_phX(p,h,X);
  //     R := if (not reducedX) then
  //       sum(data[i].R*X[i] for i in 1:size(substanceNames, 1)) else
  //       sum(data[i].R*X[i] for i in 1:size(substanceNames, 1)-1) + data[end].R*(1-sum(X[i]));
  //     d := p/(R*T);
  //   end density_phX;

     function dipoleMoment "Return mixture dynamic viscosity"
       output SI.ElectricDipoleMomentOfMolecule[nX] mu;

     algorithm
       //  for i in 1:nX loop
       // etaX[i] := fluidConstants[i].criticalTemperature;
       //  end for;

       mu := fluidConstants[:].dipoleMoment;

     end dipoleMoment;

   function critTemp "Return mixture dynamic viscosity"

    output SI.Temperature[nX] Tcrit "component dynamic viscosities";
   algorithm
   //  for i in 1:nX loop
   // etaX[i] := fluidConstants[i].criticalTemperature;
   //  end for;

   Tcrit :=fluidConstants[:].criticalTemperature;

   end critTemp;

   function critPressure "Return mixture dynamic viscosity"

    output SI.Pressure[nX] pcrit "component dynamic viscosities";
   algorithm
   //  for i in 1:nX loop
   // etaX[i] := fluidConstants[i].criticalTemperature;
   //  end for;

   pcrit :=fluidConstants[:].criticalPressure;

   end critPressure;

   function critVolume "Return mixture dynamic viscosity"

    output SI.MolarVolume[nX] Vcrit "component dynamic viscosities";
   algorithm
   //  for i in 1:nX loop
   // etaX[i] := fluidConstants[i].criticalTemperature;
   //  end for;

   Vcrit :=fluidConstants[:].criticalMolarVolume;

   end critVolume;

   function acentricFac "Return mixture dynamic viscosity"

    output Real[nX] omega "component dynamic viscosities";
   algorithm
   //  for i in 1:nX loop
   // etaX[i] := fluidConstants[i].criticalTemperature;
   //  end for;

   omega :=fluidConstants[:].acentricFactor;

   end acentricFac;

  function prandtlNumber "Return the Prandtl number"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output PrandtlNumber Pr "Prandtl number";
  algorithm
    Pr := dynamicViscosity(state)*specificHeatCapacityCp(state)/thermalConductivity(state);
  end prandtlNumber;

  function specificEnthalpy_pTX
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

  function temperature_phX "Return temperature from p, h, and X or Xi"
   extends Modelica.Icons.Function;
   input AbsolutePressure p "Pressure";
   input SpecificEnthalpy h "Specific enthalpy";
   input MassFraction X[:]=reference_X "Mass fractions";
   output Temperature T "Temperature";
  algorithm
   T := temperature(setState_phX(p,h,X));
   annotation(Documentation(info="<html></html>"));
  end temperature_phX;

  function temperature_psX "Return temperature from p,s, and X or Xi"
   extends Modelica.Icons.Function;
   input AbsolutePressure p "Pressure";
   input SpecificEntropy s "Specific entropy";
   input MassFraction X[:]=reference_X "Mass fractions";
   output Temperature T "Temperature";
  algorithm
   T := temperature(setState_psX(p,s,X));
   annotation(Documentation(info="<html></html>"));
  end temperature_psX;

  function specificEnthalpy_psX
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

  redeclare replaceable model extends FugacityCoefficient
    /*** vapour fugacity coefficient for each component ***/
    replaceable model FugacityCoeff =
        ThermalSeparation.Media.Correlations.FugacityCoefficient.IdealGas                                constrainedby ThermalSeparation.Media.Correlations.FugacityCoefficient.BaseFugacityCoefficient;
    FugacityCoeff fugacityCoeff(nS=nSubstance, T=T,p=p, y=x, v=v, Tcrit=Tcrit, pcrit=pcrit, Vcrit=Vcrit, omega=omega, mu=mu);

  equation
    phi = fugacityCoeff.phi;

  end FugacityCoefficient;
end BaseIdealGasMixture;
