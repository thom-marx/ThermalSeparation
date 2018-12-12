within ThermalSeparation.Geometry;
package SprayColumn
  import SI = Modelica.SIunits;
  record Example1
    extends Geometry(eps = 0.5, rho_solid = 1000*ones(n), c_solid = 4, a=3);

  end Example1;

  record Geometry
    extends ThermalSeparation.Geometry.BasicGeometry;

    //Packungsparameter
  //parameter Real[n] zeta = 0.01*ones(n) "pressure loss coefficient";

    //final parameter SI.Area A_free = A
      //"free cross sectional area, where vapour or liquid can be found";

  end Geometry;

  record GeometrySpray
    extends ThermalSeparation.Geometry.BasicGeometry;

    //Packungsparameter
  parameter Real[n] zeta = 0.01*ones(n) "pressure loss coefficient";

    final parameter SI.Area A_free = A
      "free cross sectional area, where vapour or liquid can be found";

  end GeometrySpray;
end SprayColumn;
