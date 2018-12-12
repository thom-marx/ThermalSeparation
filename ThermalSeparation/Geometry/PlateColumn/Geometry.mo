within ThermalSeparation.Geometry.PlateColumn;
record Geometry
  extends ThermalSeparation.Geometry.BasicGeometry;

  //Packungsparameter
  final parameter SI.VolumeFraction eps =        1
    "void fraction in the column (plate columns: ~1)";

  parameter Real h_w(max = H/n) = 0.5 "weir height";
  parameter Real l_w = 1 "weir length";
  parameter Real A_free(max = A_plate) = 0.09
    "free cross-sectional area where the vapour passes the plate";
  parameter SI.Area A_plate(max = A) = 0.9
    "area of the plate (cross-sectional area minus the area of the liquid outlet area)";
  final parameter Real phi = A_free/A_plate
    "relative freie Lochflaeche in Prozent, Stichlmair S. 8 und S. 21";
  final parameter SI.Length TS = H/n "plate spacing";
  parameter SI.Diameter d_h = 0.004 "diameter of the holes (sieve tray)";
  final parameter Real alpha = 0.61 "value of contraction (sieve tray)";
  final parameter Real zeta= (1/alpha - 1)^2 + (1-phi)^2 "s. Sti78";
 // parameter Real phi=10 "relative freie Lochflaeche in Prozent, Stichlmair S. 8 und S. 21";
end Geometry;
