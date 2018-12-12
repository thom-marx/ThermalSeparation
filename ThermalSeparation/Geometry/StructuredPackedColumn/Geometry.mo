within ThermalSeparation.Geometry.StructuredPackedColumn;
record Geometry
  extends ThermalSeparation.Geometry.BasicGeoPackedColumn;

  //Packungsparameter
 parameter Real zeta = 0.01 "pressure loss coefficient";

 //(Fluiddynamik in Füllkörper und Packungskolonnen)
  parameter Real C_L = 0.4 "constant used for pressure drop calculation";
  parameter SI.Diameter d_char = 0.01
    "characteristic length of filling material, for Re calculation";
  parameter SI.SurfaceTension sigma_crit = 33e-3
    "critical surface tension, at which the liquid is teared";
  //Onda et al.: (Mersmann)
  // PE: 33e-3 N/m, PVC: 40e-3 N/m, Keramik: 61e-3 N/m, Glas: 73e-3 N/m, Stahl: 75e-3 N/m
  parameter Real C= if d_char < 0.015 then 2 else 5.24
    "constant used for mass transfer coeff. calculation";
  //Wert aus Mersmann, p. 337
  final parameter SI.Area A_free = A*eps
    "free cross sectional area, where vapour or liquid can be found";

    //Structured Packing:
    parameter SI.Angle theta=60 "corrugation angle in deg";
    parameter SI.Length S= 0.007
    "side length of corrugation; corrugation spacing";
    final parameter Real cos_gamma = if sigma_crit < 55e-3 then 0.9 else 5.211*10^(-16.835*sigma_crit);
    final parameter Real F_SE=0.35;

    final parameter Real C1=5;
    final parameter Real C2=3;
    final parameter Real C3=0.45;
    parameter Real dp=6*(1-eps)/a;

end Geometry;
