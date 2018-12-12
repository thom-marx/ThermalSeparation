within ThermalSeparation.Geometry;
record BasicGeometry "geometry record for both column types"
    //Geometrieparameter
  parameter Integer n(min=1)=1
    "plate column: number of trays in one section, packed column: number of discrete elements in the section";
  parameter SI.SpecificHeatCapacity c_solid=100
    "specfic heat capacity of the solid in the column";
  parameter SI.Density rho_solid[n]=1000*ones(n)
    "density of the solid material in the column (NOT the column wall)";

  parameter SI.Diameter d=0.2 "inner column diameter";
 final parameter SI.Area A=Modelica.Constants.pi/4*d^2
    "Cross sectional area of the column";
  parameter SI.Length H=10 "height of the section";
  parameter SI.Length t= 0.05 "wall thickness";
  parameter Boolean hasInsulation = true "true if column is insulated";
  parameter SI.Length s= 0.01 "insulation thickness" annotation(Dialog(enable=hasInsulation));
  final parameter SI.Radius r1=d/2 "inner radius of column";
  final parameter SI.Radius r2=r1+t
    "outer radius of column (not considering insulation)";
  final parameter SI.Radius r3=if hasInsulation then r1+t+s else r1+t
    "outer radius of column + insulation";
  parameter SI.Density rho_wall = 1000 "density of the column wall";
  parameter SI.SpecificHeatCapacity c_wall = 1000
    "specific heat capacity of the column wall";
  parameter SI.ThermalConductivity lambda_wall = 50
    "Thermal conductivity of the column wall";
  parameter SI.ThermalConductivity lambda_ins = 0.02
    "Thermal conductivity of insulation" annotation(Dialog(enable=hasInsulation));
end BasicGeometry;
