within ThermalSeparation.Media;
package BaseMediumVapour "base class for vapour medium package"
  constant Integer nSubstance(min=2)=2 "number of substances";
  constant SI.Volume V[nSubstance] "diffusion volumes, Taylor p. 69";
  constant Integer eq_Tsonopoulos[nSubstance]
    "no. of equation for equation a and b in Table 4-5 in Properties of Gases and Liquids, 5th ed.";
    constant Real sigma[nSubstance] "surface tension";
    constant Real epsilon_k[nSubstance]
    "Lennard-Jones force constant (epsilon/k), constant value";
    constant SI.MolarMass[nSubstance] MMX "molar masses of components";
    constant Boolean delta_hv_medium = true
    "true, if evaporation enthalpy is considered in medium model";

    record ThermodynamicState "thermodynamic state variables"
    extends Modelica.Icons.Record;
    ThermalSeparation.Media.Types.AbsolutePressure p
      "Absolute pressure of medium";
    ThermalSeparation.Media.Types.Temperature T "Temperature of medium";
    ThermalSeparation.Media.Types.MassFraction X[
                   nSubstance]
      "Mass fractions (= (component mass)/total mass  m_i/m)";
    annotation(Documentation(info="<html></html>"));
    end ThermodynamicState;

  replaceable partial model BaseProperties
    "Base properties (p, d, T, h, u, R, MM and X) of a medium"

    input ThermalSeparation.Media.Types.AbsolutePressure p;
    input SI.MoleFraction x[nSubstance];
    input ThermalSeparation.Media.Types.Temperature T;
    input SI.Concentration c[nSubstance];
    input SI.MoleFraction x_star[nSubstance];

    parameter SI.Temperature T0 = 293.15 "reference temperature";

    input ThermalSeparation.Units.MolarEnthalpy h
      "Specific enthalpy of medium";
    output SI.MolarInternalEnergy u "Specific internal energy of medium";
    output ThermalSeparation.Media.Types.Density d "mixture density";
    output ThermalSeparation.Media.Types.MolarMass MM "mixture molar mass";
    output SI.DynamicViscosity eta;
    output SI.ThermalConductivity lambda;
    output SI.SpecificHeatCapacity cp;

    output ThermodynamicState state "thermodynamic state variables: p, T, X";

    output ThermodynamicProperties properties;

    SI.MassFraction[nSubstance] X
      "Mass fractions (= (component mass)/total mass  m_i/m)";

    /*** for feed stage ***/
    ThermalSeparation.Interfaces.MediumConOut mediumConOut(h=h, rho = d, MM=MM);

  equation
    // connect state with BaseProperties
    state.T = T;
    state.p = p;
    state.X = X;

    properties.T=T;
    properties.eta=eta;
    properties.rho=d;
    properties.MM = MM;
    properties.x = x;
    properties.lambda = lambda;
    properties.cp = cp;
    properties.u = u;
    properties.c = c;
    properties.d=d;
    properties.h=h;

    annotation (structurallyIncomplete,
           Icon(graphics={Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid), Text(extent={{-152,164},{152,
                  102}}, textString =    "%name")}));
  end BaseProperties;

  replaceable partial model CalcSpecificEnthalpy
  /*** Enthalpy ***/
    input SI.Temperature T;
    input SI.MoleFraction x[nSubstance];
    input SI.Pressure p;
    parameter SI.Temperature T0 = 293.15 "reference temperature";

    output ThermalSeparation.Units.MolarEnthalpy h;

    annotation (Icon(graphics={
                          Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end CalcSpecificEnthalpy;

  replaceable partial model EvaporationEnthalpy
    /*** evaporation enthalpy for each component ***/
    input SI.Pressure p;
    input SI.Temperature T;
    output ThermalSeparation.Units.MolarEnthalpy h[nSubstance];

  annotation (Icon(graphics={
                        Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end EvaporationEnthalpy;

  replaceable partial model FugacityCoefficient
    /*** vapour fugacity coefficient for each component ***/

    input SI.Pressure p;
    input SI.Temperature T;
    input SI.MoleFraction x[nSubstance];
    input SI.MolarVolume v;
    output Real phi[nSubstance];

  annotation (Icon(graphics={
                        Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end FugacityCoefficient;

replaceable model DiffusionCoefficient
    "returns the binary diffusion coefficients of the components in the mixture"

   input SI.Temperature T;
   input SI.Pressure p;

   // Diffusion coefficents in the following form: D12, D13, D14, D23, D24, D34
   // the numbers correspond to the ordering of the substances
   output SI.DiffusionCoefficient D[a] "Maxwell-Stefan diffusion coefficients";
   output SI.DiffusionCoefficient D_matrix[nSubstance,nSubstance]
      "Maxwell-Stefan diffusion coefficients in matrix form, where D_12 = D_21";
   output SI.DiffusionCoefficient D_diluted[nSubstance]
      "diffusion coefficient of a component in the mixture";

  protected
  parameter Integer aux[:] = {1,3,6,10,15, 21, 28, 36, 45};
  parameter Integer a = aux[nSubstance-1]
      "number of binary diffusion coefficients depending on the number of substances";
replaceable model DiffCoeff =
    ThermalSeparation.Media.Correlations.DiffusionCoefficient.Vapour.Fuller;

/*** model where the diffusion coefficients in a binary mixture are used to calculate the diffusion coefficients in a multicomponent mixture ***/
DiffCoeff diffCoeff(T=T, p=p, nS=nSubstance, MMX=MMX, V=V);

equation
D = diffCoeff.D;
for i in 1:nSubstance loop
D_diluted[i] = sum(diffCoeff.D_matrix[i,:])/nSubstance;
end for;

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

  replaceable model ThermodynamicFactor

    output Real Gamma[nSubstance,nSubstance] = diagonal(ones(nSubstance));
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

  SI.DynamicViscosity eta;

  ThermalSeparation.Media.Types.Density rho "Density of medium";

  ThermalSeparation.Media.Types.MolarMass MM
      "Molar mass (of mixture or single fluid)";

  SI.MoleFraction x[nSubstance];
  SI.Density d;
  ThermalSeparation.Units.MolarEnthalpy h;

  SI.ThermalConductivity lambda;
  SI.SpecificHeatCapacity cp;
  SI.MolarInternalEnergy u;
  SI.Concentration c[nSubstance];

end ThermodynamicProperties;
end BaseMediumVapour;
