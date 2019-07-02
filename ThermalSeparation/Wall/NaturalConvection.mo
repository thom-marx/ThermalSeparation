within ThermalSeparation.Wall;
model NaturalConvection
  "Heat Transfer with natural convection of the outside air"
  extends ThermalSeparation.Wall.BaseWall;
  parameter SI.Pressure p_amb = 1e5 "ambient pressure";
   SI.SurfaceCoefficientOfHeatTransfer alpha[n](start=0.1*ones(n))
    "Surface Coefficient of Heat Transfer at outer Wall of the column";
protected
  SI.NusseltNumber Nu[n] "Nusselt Number";
  SI.RayleighNumber Ra[n] "Rayleigh Number";
  SI.GrashofNumber Gr[n] "Grasshof Number";
  SI.PrandtlNumber Pr "Prandtl Number";
  package Medium = Modelica.Media.Air.SimpleAir
    "Simple Air to calculate the dimensionless numbers at ambient Temperature and Pressure";
  Medium.BaseProperties air(p=p_amb,T=T_amb);
  Medium.SpecificHeatCapacity cp=Medium.specificHeatCapacityCp(air);
  Medium.DynamicViscosity eta=Medium.dynamicViscosity(air);
  Medium.Density rho=Medium.density(air);
  Medium.ThermalConductivity lambda_air = Medium.thermalConductivity(air);
algorithm
  Pr :=eta*cp/lambda_air;
  for j in 1:n loop
    Gr[j] :=Modelica.Constants.g_n*geometry.H^3*rho/eta*abs(T_wall[j] - T_amb)/T_amb;
    Ra[j] :=Gr[j]*Pr;
    Nu[j] :=(0.825 + 0.387*(Ra[j]*(1 + (0.492/Pr)^(9/16))^(-16/9))^(1/6))^2 + 0.435*
    geometry.H/(2*geometry.r2);
    alpha[j] :=Nu[j]*lambda_air/geometry.H;
    if geometry.hasInsulation then
      Qdot_out[j]:=1/(1/alpha[j] + geometry.s/geometry.lambda_ins)*2*Modelica.Constants.pi*geometry.H/n*(geometry.r2 + geometry.r3)/2*(T_wall[j] - T_amb);
    else
      Qdot_out[j]:=alpha[j]*2*Modelica.Constants.pi*geometry.H/n*(geometry.r2 + geometry.r3)/2*(T_wall[j] - T_amb);
    end if;
  end for;
  annotation (Icon(graphics),
                    DymolaStoredErrors,
    Documentation(info="<html> This model calculates the heat loss to ambiance using a correlation for the outside heat transfer coefficient of air obtained from 
VDI Wärmeatlas, 6.Auflage, 1991, S. Fa2. For further information of this model see also  <a href=\"Modelica://ThermalSeparation.Wall.BaseWall\">BaseWall</a></html>"));
end NaturalConvection;
