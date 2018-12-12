within ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common;
partial package SingleGasNasa
  "Medium model of an ideal gas based on NASA source"

  extends
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialPureSubstance(
    mediumName=data.name,
    substanceNames={data.name},
    final reducedX=true,
    singleState=false,
    Temperature(
      min=200,
      max=6000,
      start=500,
      nominal=500),
    SpecificEnthalpy(start=if referenceChoice == ReferenceEnthalpy.ZeroAt0K then
                data.H0 else if referenceChoice == ReferenceEnthalpy.UserDefined then
                h_offset else 0, nominal=1.0e5),
    Density(start=10, nominal=10),
    AbsolutePressure(start=10e5, nominal=10e5));

  redeclare replaceable record extends ThermodynamicState
    "thermodynamic state variables for ideal gases"
    AbsolutePressure p "Absolute pressure of medium";
    Temperature T "Temperature of medium";
  end ThermodynamicState;

  redeclare replaceable record extends FluidConstants
    "Extended fluid constants"
    Temperature criticalTemperature "critical temperature";
    AbsolutePressure criticalPressure "critical pressure";
    MolarVolume criticalMolarVolume "critical molar Volume";
    Real acentricFactor "Pitzer acentric factor";
    Temperature triplePointTemperature "triple point temperature";
    AbsolutePressure triplePointPressure "triple point pressure";
    Temperature meltingPoint "melting point at 101325 Pa";
    Temperature normalBoilingPoint "normal boiling point (at 101325 Pa)";
    Modelica.SIunits.ElectricDipoleMomentOfMolecule dipoleMoment
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

  import SI = Modelica.SIunits;
  import Modelica.Math;
  import
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices.ReferenceEnthalpy;

  constant Boolean excludeEnthalpyOfFormation=true
    "If true, enthalpy of formation Hf is not included in specific enthalpy h";
  constant ReferenceEnthalpy.Temp referenceChoice=Choices.
        ReferenceEnthalpy.ZeroAt0K "Choice of reference enthalpy";
  constant SpecificEnthalpy h_offset=0.0
    "User defined offset for reference enthalpy, if referenceChoice = UserDefined";

  constant
    ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord data
    "Data record of ideal gas substance";

  constant FluidConstants[nS] fluidConstants "constant data for the fluid";

  redeclare model extends BaseProperties(
   T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
   p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default))
    "Base properties of ideal gas medium"

  equation
    //     assert(T >= 200 and T <= 6000, "
    // Temperature T (= " + String(T) + " K) is not in the allowed range
    // 200 K <= T <= 6000 K required from medium model \"" + mediumName + "\".
    // ");
    MM = data.MM;
    R = data.R;
    h = h_T(data, T, excludeEnthalpyOfFormation, referenceChoice, h_offset);
    u = h - R*T;

    // Has to be written in the form d=f(p,T) in order that static
    // state selection for p and T is possible
    d = p/(R*T);
    // connect state with BaseProperties
    state.T = T;
    state.p = p;
   annotation(structurallyIncomplete);
  end BaseProperties;

    redeclare function setState_pTX
    "Return thermodynamic state as function of p, T and composition X"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input Temperature T "Temperature";
      input MassFraction X[:]=reference_X "Mass fractions";
      output ThermodynamicState state;
    algorithm
      state := ThermodynamicState(p=p,T=T);
    end setState_pTX;

    redeclare function setState_phX
    "Return thermodynamic state as function of p, h and composition X"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input SpecificEnthalpy h "Specific enthalpy";
      input MassFraction X[:]=reference_X "Mass fractions";
      output ThermodynamicState state;
    algorithm
      state := ThermodynamicState(p=p,T=T_h(h));
    end setState_phX;

    redeclare function setState_psX
    "Return thermodynamic state as function of p, s and composition X"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input SpecificEntropy s "Specific entropy";
      input MassFraction X[:]=reference_X "Mass fractions";
      output ThermodynamicState state;
    algorithm
      state := ThermodynamicState(p=p,T=T_ps(p,s));
    end setState_psX;

    redeclare function setState_dTX
    "Return thermodynamic state as function of d, T and composition X"
      extends Modelica.Icons.Function;
      input Density d "density";
      input Temperature T "Temperature";
      input MassFraction X[:]=reference_X "Mass fractions";
      output ThermodynamicState state;
    algorithm
      state := ThermodynamicState(p=d*data.R*T,T=T);
    end setState_dTX;

  redeclare function extends pressure "return pressure of ideal gas"
  algorithm
    p := state.p;
  end pressure;

  redeclare function extends temperature "return temperature of ideal gas"
  algorithm
    T := state.T;
  end temperature;

  redeclare function extends density "return density of ideal gas"
  algorithm
    d := state.p/(data.R*state.T);
  end density;

  redeclare function extends specificEnthalpy "Return specific enthalpy"
    extends Modelica.Icons.Function;
  algorithm
    h := h_T(data,state.T);
  end specificEnthalpy;

  redeclare function extends specificInternalEnergy
    "Return specific internal energy"
    extends Modelica.Icons.Function;
  algorithm
    u := h_T(data,state.T) - data.R*state.T;
  end specificInternalEnergy;

  redeclare function extends specificEntropy "Return specific entropy"
    extends Modelica.Icons.Function;
  algorithm
    s := s0_T(data, state.T) - data.R*Modelica.Math.log(state.p/reference_p);
  end specificEntropy;

  redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
    extends Modelica.Icons.Function;
  algorithm
    g := h_T(data,state.T) - state.T*specificEntropy(state);
  end specificGibbsEnergy;

  redeclare function extends specificHelmholtzEnergy
    "Return specific Helmholtz energy"
    extends Modelica.Icons.Function;
  algorithm
    f := h_T(data,state.T) - data.R*state.T - state.T*specificEntropy(state);
  end specificHelmholtzEnergy;

  redeclare function extends specificHeatCapacityCp
    "Return specific heat capacity at constant pressure"
  algorithm
    cp := cp_T(data, state.T);
  end specificHeatCapacityCp;

  redeclare function extends specificHeatCapacityCv
    "Compute specific heat capacity at constant volume from temperature and gas data"
  algorithm
    cv := cp_T(data, state.T) - data.R;
  end specificHeatCapacityCv;

  redeclare function extends isentropicExponent "Return isentropic exponent"
  algorithm
    gamma := specificHeatCapacityCp(state)/specificHeatCapacityCv(state);
  end isentropicExponent;

  redeclare function extends velocityOfSound "Return velocity of sound"
    extends Modelica.Icons.Function;
  algorithm
    a := sqrt(data.R*state.T*cp_T(data, state.T)/specificHeatCapacityCv(state));
  end velocityOfSound;

  function isentropicEnthalpyApproximation
    "approximate method of calculating h_is from upstream properties and downstream pressure"
    extends Modelica.Icons.Function;
    input SI.Pressure p2 "downstream pressure";
    input ThermodynamicState state "properties at upstream location";
    input Boolean exclEnthForm=excludeEnthalpyOfFormation
      "If true, enthalpy of formation Hf is not included in specific enthalpy h";
    input ReferenceEnthalpy.Temp refChoice=referenceChoice
      "Choice of reference enthalpy";
    input SpecificEnthalpy h_off=h_offset
      "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
    output SI.SpecificEnthalpy h_is "isentropic enthalpy";
  protected
    IsentropicExponent gamma =  isentropicExponent(state) "Isentropic exponent";
  algorithm
    h_is := h_T(data,state.T,exclEnthForm,refChoice,h_off) +
      gamma/(gamma - 1.0)*state.p/density(state)*((p2/state.p)^((gamma - 1)/gamma) - 1.0);
  end isentropicEnthalpyApproximation;

  redeclare function extends isentropicEnthalpy "Return isentropic enthalpy"
  input Boolean exclEnthForm=excludeEnthalpyOfFormation
      "If true, enthalpy of formation Hf is not included in specific enthalpy h";
  input ReferenceEnthalpy.Temp refChoice=referenceChoice
      "Choice of reference enthalpy";
  input SpecificEnthalpy h_off=h_offset
      "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
  algorithm
    h_is := isentropicEnthalpyApproximation(p_downstream,refState,exclEnthForm,refChoice,h_off);
  end isentropicEnthalpy;

  redeclare function extends isobaricExpansionCoefficient
    "Returns overall the isobaric expansion coefficient beta"
  algorithm
    beta := 1/state.T;
  end isobaricExpansionCoefficient;

  redeclare function extends isothermalCompressibility
    "Returns overall the isothermal compressibility factor"
  algorithm
    kappa := 1.0/state.p;
  end isothermalCompressibility;

  redeclare function extends density_derp_T
    "Returns the partial derivative of density with respect to pressure at constant temperature"
  algorithm
    ddpT := 1/(state.T*data.R);
  end density_derp_T;

  redeclare function extends density_derT_p
    "Returns the partial derivative of density with respect to temperature at constant pressure"
  algorithm
    ddTp := -state.p/(state.T*state.T*data.R);
  end density_derT_p;

  redeclare function extends density_derX
    "Returns the partial derivative of density with respect to mass fractions at constant pressure and temperature"
  algorithm
    dddX := fill(0,0);
  end density_derX;

  function cp_T
    "Compute specific heat capacity at constant pressure from temperature and gas data"
    extends Modelica.Icons.Function;
    input
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord
      data "Ideal gas data";
    input SI.Temperature T "Temperature";
    output SI.SpecificHeatCapacity cp "Specific heat capacity at temperature T";
  algorithm
    cp := smooth(0,if T < data.Tlimit then data.R*(1/(T*T)*(data.alow[1] + T*(
      data.alow[2] + T*(1.*data.alow[3] + T*(data.alow[4] + T*(data.alow[5] + T
      *(data.alow[6] + data.alow[7]*T))))))) else data.R*(1/(T*T)*(data.ahigh[1]
       + T*(data.ahigh[2] + T*(1.*data.ahigh[3] + T*(data.ahigh[4] + T*(data.
      ahigh[5] + T*(data.ahigh[6] + data.ahigh[7]*T))))))));
    annotation (InlineNoEvent=false);
  end cp_T;

  function cp_Tlow
    "Compute specific heat capacity at constant pressure, low T region"
    extends Modelica.Icons.Function;
    input
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord
      data "Ideal gas data";
    input SI.Temperature T "Temperature";
    output SI.SpecificHeatCapacity cp "Specific heat capacity at temperature T";
  algorithm
    cp := data.R*(1/(T*T)*(data.alow[1] + T*(
      data.alow[2] + T*(1.*data.alow[3] + T*(data.alow[4] + T*(data.alow[5] + T
      *(data.alow[6] + data.alow[7]*T)))))));
    annotation (Inline=false, derivative(zeroDerivative=data) = cp_Tlow_der);
  end cp_Tlow;

  function cp_Tlow_der
    "Compute specific heat capacity at constant pressure, low T region"
    extends Modelica.Icons.Function;
    input
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord
      data "Ideal gas data";
    input SI.Temperature T "Temperature";
    input Real dT "Temperature derivative";
    output Real cp_der "Derivative of specific heat capacity";
  algorithm
    cp_der := dT*data.R/(T*T*T)*(-2*data.alow[1] + T*(
      -data.alow[2] + T*T*(data.alow[4] + T*(2.*data.alow[5] + T
      *(3.*data.alow[6] + 4.*data.alow[7]*T)))));
  end cp_Tlow_der;

  function h_T "Compute specific enthalpy from temperature and gas data; reference is decided by the 
    refChoice input, or by the referenceChoice package constant by default"
    import
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices;
    extends Modelica.Icons.Function;
    input
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord
      data "Ideal gas data";
    input SI.Temperature T "Temperature";
    input Boolean exclEnthForm=excludeEnthalpyOfFormation
      "If true, enthalpy of formation Hf is not included in specific enthalpy h";
    input Choices.ReferenceEnthalpy.Temp refChoice =                                                              referenceChoice
      "Choice of reference enthalpy";
    input SI.SpecificEnthalpy h_off=h_offset
      "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
    output SI.SpecificEnthalpy h "Specific enthalpy at temperature T";
      //     annotation (InlineNoEvent=false, Inline=false,
      //                 derivative(zeroDerivative=data,
      //                            zeroDerivative=exclEnthForm,
      //                            zeroDerivative=refChoice,
      //                            zeroDerivative=h_off) = h_T_der);
  algorithm
    h := smooth(0,(if T < data.Tlimit then data.R*((-data.alow[1] + T*(data.
      blow[1] + data.alow[2]*Math.log(T) + T*(1.*data.alow[3] + T*(0.5*data.
      alow[4] + T*(1/3*data.alow[5] + T*(0.25*data.alow[6] + 0.2*data.alow[7]*T))))))
      /T) else data.R*((-data.ahigh[1] + T*(data.bhigh[1] + data.ahigh[2]*
      Math.log(T) + T*(1.*data.ahigh[3] + T*(0.5*data.ahigh[4] + T*(1/3*data.
      ahigh[5] + T*(0.25*data.ahigh[6] + 0.2*data.ahigh[7]*T))))))/T)) + (if
      exclEnthForm then -data.Hf else 0.0) + (if (refChoice
       == Choices.ReferenceEnthalpy.ZeroAt0K) then data.H0 else 0.0) + (if
      refChoice == Choices.ReferenceEnthalpy.UserDefined then h_off else
            0.0));

    annotation (Inline=false,smoothOrder=1);
  end h_T;

  function h_T_der "derivative function for h_T"
    import
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices;
    extends Modelica.Icons.Function;
    input
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord
      data "Ideal gas data";
    input SI.Temperature T "Temperature";
    input Boolean exclEnthForm=excludeEnthalpyOfFormation
      "If true, enthalpy of formation Hf is not included in specific enthalpy h";
    input Choices.ReferenceEnthalpy.Temp refChoice=referenceChoice
      "Choice of reference enthalpy";
    input SI.SpecificEnthalpy h_off=h_offset
      "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
    input Real dT "Temperature derivative";
    output Real h_der "Specific enthalpy at temperature T";
  algorithm
    h_der := dT*cp_T(data,T);
  end h_T_der;

  function h_Tlow "Compute specific enthalpy, low T region; reference is decided by the 
    refChoice input, or by the referenceChoice package constant by default"
    import
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices;
    extends Modelica.Icons.Function;
    input
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord
      data "Ideal gas data";
    input SI.Temperature T "Temperature";
    input Boolean exclEnthForm=excludeEnthalpyOfFormation
      "If true, enthalpy of formation Hf is not included in specific enthalpy h";
    input Choices.ReferenceEnthalpy.Temp refChoice=referenceChoice
      "Choice of reference enthalpy";
    input SI.SpecificEnthalpy h_off=h_offset
      "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
    output SI.SpecificEnthalpy h "Specific enthalpy at temperature T";
      //     annotation (Inline=false,InlineNoEvent=false, derivative(zeroDerivative=data,
      //                                zeroDerivative=exclEnthForm,
      //                                zeroDerivative=refChoice,
      //                                zeroDerivative=h_off) = h_Tlow_der);
  algorithm
    h := data.R*((-data.alow[1] + T*(data.
      blow[1] + data.alow[2]*Math.log(T) + T*(1.*data.alow[3] + T*(0.5*data.
      alow[4] + T*(1/3*data.alow[5] + T*(0.25*data.alow[6] + 0.2*data.alow[7]*T))))))
      /T) + (if
      exclEnthForm then -data.Hf else 0.0) + (if (refChoice
       == Choices.ReferenceEnthalpy.ZeroAt0K) then data.H0 else 0.0) + (if
      refChoice == Choices.ReferenceEnthalpy.UserDefined then h_off else
            0.0);
    annotation(Inline=false,InlineNoEvent=false,smoothOrder=1);
  end h_Tlow;

  function h_Tlow_der "Compute specific enthalpy, low T region; reference is decided by the 
    refChoice input, or by the referenceChoice package constant by default"
    import
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.PartialMedium.Choices;
    extends Modelica.Icons.Function;
    input
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord
      data "Ideal gas data";
    input SI.Temperature T "Temperature";
    input Boolean exclEnthForm=excludeEnthalpyOfFormation
      "If true, enthalpy of formation Hf is not included in specific enthalpy h";
    input Choices.ReferenceEnthalpy.Temp refChoice=referenceChoice
      "Choice of reference enthalpy";
    input SI.SpecificEnthalpy h_off=h_offset
      "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
    input Real dT(unit="K/s") "Temperature derivative";
    output Real h_der(unit="J/(kg.s)")
      "Derivative of specific enthalpy at temperature T";
  algorithm
    h_der := dT*cp_Tlow(data,T);
  end h_Tlow_der;

  function s0_T "Compute specific entropy from temperature and gas data"
    extends Modelica.Icons.Function;
    input
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord
      data "Ideal gas data";
    input SI.Temperature T "Temperature";
    output SI.SpecificEntropy s "Specific entropy at temperature T";
    //    annotation (InlineNoEvent=false);
  algorithm
    s := noEvent(if T < data.Tlimit then data.R*(data.blow[2] - 0.5*data.alow[
      1]/(T*T) - data.alow[2]/T + data.alow[3]*Math.log(T) + T*(
      data.alow[4] + T*(0.5*data.alow[5] + T*(1/3*data.alow[6] + 0.25*data.alow[
      7]*T)))) else data.R*(data.bhigh[2] - 0.5*data.ahigh[1]/(T*T) - data.
      ahigh[2]/T + data.ahigh[3]*Math.log(T) + T*(data.ahigh[4]
       + T*(0.5*data.ahigh[5] + T*(1/3*data.ahigh[6] + 0.25*data.ahigh[7]*T)))));
  end s0_T;

  function s0_Tlow "Compute specific entropy, low T region"
    extends Modelica.Icons.Function;
    input
      ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord
      data "Ideal gas data";
    input SI.Temperature T "Temperature";
    output SI.SpecificEntropy s "Specific entropy at temperature T";
    //    annotation (InlineNoEvent=false);
  algorithm
    s := data.R*(data.blow[2] - 0.5*data.alow[
      1]/(T*T) - data.alow[2]/T + data.alow[3]*Math.log(T) + T*(
      data.alow[4] + T*(0.5*data.alow[5] + T*(1/3*data.alow[6] + 0.25*data.alow[
      7]*T))));
  end s0_Tlow;

  function dynamicViscosityLowPressure
    "Dynamic viscosity of low pressure gases"
    extends Modelica.Icons.Function;
    input SI.Temp_K T "Gas temperature";
    input SI.Temp_K Tc "Critical temperature of gas";
    input SI.MolarMass M "Molar mass of gas";
    input SI.MolarVolume Vc "Critical molar volume of gas";
    input Real w "Acentric factor of gas";
    input Real mu "Dipole moment of gas molecule";
    input Real k =  0.0 "Special correction for highly polar substances";
    output SI.DynamicViscosity eta "Dynamic viscosity of gas";
  protected
    parameter Real Const1_SI=40.785*10^(-9.5)
      "Constant in formula for eta converted to SI units";
    parameter Real Const2_SI=131.3/1000.0
      "Constant in formula for mur converted to SI units";
    parameter Real mur=Const2_SI*mu/sqrt(Vc*Tc)
      "Dimensionless dipole moment of gas molecule";
    parameter Real Fc=1 - 0.2756*w + 0.059035*mur^4 + k
      "Factor to account for molecular shape and polarities of gas";
    Real Tstar "Dimensionless temperature defined by equation below";
    Real Ov "Viscosity collision integral for the gas";

  algorithm
    Tstar := 1.2593*T/Tc;
    Ov := 1.16145*Tstar^(-0.14874) + 0.52487*exp(-0.7732*Tstar) + 2.16178*exp(-2.43787
      *Tstar);
    eta := Const1_SI*Fc*sqrt(M*T)/(Vc^(2/3)*Ov);
    annotation (Documentation(info="<html>
<p>
The used formula are based on the method of Chung et al (1984, 1988) referred to in ref [1] chapter 9.
The formula 9-4.10 is the one being used. The Formula is given in non-SI units, the follwong onversion constants were used to
transform the formula to SI units:
</p>
 
<ul>
<li> <b>Const1_SI:</b> The factor 10^(-9.5) =10^(-2.5)*1e-7 where the 
     factor 10^(-2.5) originates from the conversion of g/mol->kg/mol + cm^3/mol->m^3/mol
      and the factor 1e-7 is due to conversionfrom microPoise->Pa.s.</li>
<li>  <b>Const2_SI:</b> The factor 1/3.335641e-27 = 1e-3/3.335641e-30 
      where the factor 3.335641e-30 comes from debye->C.m and
      1e-3 is due to conversion from cm^3/mol->m^3/mol</li>
</ul>
 
<h4>References:</h4>
<p>
[1] Bruce E. Poling, John E. Prausnitz, John P. O'Connell, \"The Properties of Gases and Liquids\" 5th Ed. Mc Graw Hill.
</p>
 
<h4>Author</h4>
<p>T. Skoglund, Lund, Sweden, 2004-08-31</p>
 
</html>"));
  end dynamicViscosityLowPressure;

  redeclare replaceable function extends dynamicViscosity "dynamic viscosity"
  algorithm
    assert(fluidConstants[1].hasCriticalData,
    "Failed to compute dynamicViscosity: For the species \"" + mediumName + "\" no critical data is available.");
    assert(fluidConstants[1].hasDipoleMoment,
    "Failed to compute dynamicViscosity: For the species \"" + mediumName + "\" no critical data is available.");
    eta := dynamicViscosityLowPressure(state.T,
                       fluidConstants[1].criticalTemperature,
                       fluidConstants[1].molarMass,
                       fluidConstants[1].criticalMolarVolume,
                       fluidConstants[1].acentricFactor,
                       fluidConstants[1].dipoleMoment);
  end dynamicViscosity;

  function thermalConductivityEstimate
    "Thermal conductivity of polyatomic gases(Eucken and Modified Eucken correlation)"
    extends Modelica.Icons.Function;
    input SpecificHeatCapacity Cp "Constant pressure heat capacity";
    input DynamicViscosity eta "Dynamic viscosity";
    input Integer method(min=1,max=2)=1
      "1: Eucken Method, 2: Modified Eucken Method";
    output ThermalConductivity lambda "Thermal conductivity [W/(m.k)]";
  algorithm
    lambda := if method == 1 then eta*(Cp - data.R + (9/4)*data.R) else eta*(Cp
       - data.R)*(1.32 + 1.77/((Cp/Modelica.Constants.R) - 1.0));
    annotation (Documentation(info="<html>
<p>
This function provides two similar methods for estimating the 
thermal conductivity of polyatomic gases.
The Eucken method (input method == 1) gives good results for low temperatures, 
but it tends to give an underestimated value of the thermal conductivity 
(lambda) at higher temperatures.<br>
The Modified Eucken method (input method == 2) gives good results for 
high-temperatures, but it tends to give an overestimated value of the
thermal conductivity (lambda) at low temperatures.
</p>
</html>"));
  end thermalConductivityEstimate;

  redeclare replaceable function extends thermalConductivity
    "thermal conductivity of gas"
    input Integer method=1 "1: Eucken Method, 2: Modified Eucken Method";
  algorithm
    assert(fluidConstants[1].hasCriticalData,
    "Failed to compute thermalConductivity: For the species \"" + mediumName + "\" no critical data is available.");
    lambda := thermalConductivityEstimate(specificHeatCapacityCp(state),
      dynamicViscosity(state), method=method);
  end thermalConductivity;

  redeclare function extends molarMass "return the molar mass of the medium"
  algorithm
    MM := data.MM;
  end molarMass;

  function T_h "Compute temperature from specific enthalpy"
    input SpecificEnthalpy h "Specific enthalpy";
    output Temperature T "Temperature";
  protected
  package Internal
      "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
    extends
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.OneNonLinearEquation;
    redeclare record extends f_nonlinear_Data
        "Data to be passed to non-linear function"
      extends
          ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord;
    end f_nonlinear_Data;

    redeclare function extends f_nonlinear
    algorithm
        y := h_T(f_nonlinear_data,x);
    end f_nonlinear;

    // Dummy definition has to be added for current Dymola
    redeclare function extends solve
    end solve;
  end Internal;

  algorithm
    T := Internal.solve(h, 200, 6000, 1.0e5, {1}, data);
  end T_h;

  function T_ps "Compute temperature from pressure and specific entropy"
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    output Temperature T "Temperature";
  protected
  package Internal
      "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
    extends
        ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.OneNonLinearEquation;
    redeclare record extends f_nonlinear_Data
        "Data to be passed to non-linear function"
      extends
          ThermalSeparation.Media.IdealGasMixtures.BaseClasses.Common.DataRecord;
    end f_nonlinear_Data;

    redeclare function extends f_nonlinear
    algorithm
        y := s0_T(f_nonlinear_data,x)- data.R*Modelica.Math.log(p/reference_p);
    end f_nonlinear;

    // Dummy definition has to be added for current Dymola
    redeclare function extends solve
    end solve;
  end Internal;

  algorithm
    T := Internal.solve(s, 200, 6000, p, {1}, data);
  end T_ps;

  annotation (
    Documentation(info="<HTML>
<p>
This model calculates medium properties
for an ideal gas of a single substance, or for an ideal
gas consisting of several substances where the
mass fractions are fixed. Independent variables
are temperature <b>T</b> and pressure <b>p</b>.
Only density is a function of T and p. All other quantities
are solely a function of T. The properties
are valid in the range:
</p>
<pre>
   200 K &le; T &le; 6000 K
</pre>
<p>
The following quantities are always computed:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td valign=\"top\"><b>Variable</b></td>
      <td valign=\"top\"><b>Unit</b></td>
      <td valign=\"top\"><b>Description</b></td></tr>
  <tr><td valign=\"top\">h</td>
      <td valign=\"top\">J/kg</td>
      <td valign=\"top\">specific enthalpy h = h(T)</td></tr>
  <tr><td valign=\"top\">u</td>
      <td valign=\"top\">J/kg</td>
      <td valign=\"top\">specific internal energy u = u(T)</b></td></tr>
  <tr><td valign=\"top\">d</td>
      <td valign=\"top\">kg/m^3</td>
      <td valign=\"top\">density d = d(p,T)</td></tr>
</table>
<p>
For the other variables, see the functions in
Modelica.Media.IdealGases.Common.SingleGasNasa.
Note, dynamic viscosity and thermal conductivity are only provided
for gases that use a data record from Modelica.Media.IdealGases.FluidData.
Currently these are the following gases:
</p>
<pre>
  Ar
  C2H2_vinylidene
  C2H4
  C2H5OH
  C2H6
  C3H6_propylene
  C3H7OH
  C3H8
  C4H8_1_butene
  C4H9OH
  C4H10_n_butane
  C5H10_1_pentene
  C5H12_n_pentane
  C6H6
  C6H12_1_hexene
  C6H14_n_heptane
  C7H14_1_heptene
  C8H10_ethylbenz
  CH3OH
  CH4
  CL2
  CO
  CO2
  F2
  H2
  H2O
  He
  N2
  N2O 
  NH3
  NO
  O2
  SO2
  SO3
</pre>
<p>
<b>Sources for model and literature:</b><br>
Original Data: Computer program for calculation of complex chemical
equilibrium compositions and applications. Part 1: Analysis
Document ID: 19950013764 N (95N20180) File Series: NASA Technical Reports
Report Number: NASA-RP-1311  E-8017  NAS 1.61:1311
Authors: Gordon, Sanford (NASA Lewis Research Center)
 Mcbride, Bonnie J. (NASA Lewis Research Center)
Published: Oct 01, 1994.
</p>
<p><b>Known limits of validity:</b></br>
The data is valid for
temperatures between 200K and 6000K.  A few of the data sets for
monatomic gases have a discontinuous 1st derivative at 1000K, but
this never caused problems so far.
</p>
<p>
This model has been copied from the ThermoFluid library
and adapted to the Modelica.Media package.
</p>
</HTML>"),
    Icon(graphics),
    Diagram(graphics));
end SingleGasNasa;
