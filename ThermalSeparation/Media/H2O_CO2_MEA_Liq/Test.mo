within ThermalSeparation.Media.H2O_CO2_MEA_Liq;
model Test

parameter Integer nS=3;

parameter Modelica.Units.SI.Temperature T=383;
parameter Modelica.Units.SI.Temperature T0=293;
parameter Modelica.Units.SI.Pressure p=2e5;
parameter Modelica.Units.SI.MoleFraction x[nS]={0.888,0.037,0.075};

/* dummy inouts */

BaseProperties baseProperties(T=T, p=p, x=x);

CalcSpecificEnthalpy calcSpecificEnthalpy(T0=T0, T=T, p=p, x=x);

end Test;
