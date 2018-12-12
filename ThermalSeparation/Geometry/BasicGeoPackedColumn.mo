within ThermalSeparation.Geometry;
record BasicGeoPackedColumn "base record for packed column"
  extends BasicGeometry;
  parameter SI.VolumeFraction eps =         0.5 "void fraction in the column";
   parameter ThermalSeparation.Units.VolumetricArea a=
                                         600
    "specific gas-liquid interfacial area in m2/m3";

end BasicGeoPackedColumn;
