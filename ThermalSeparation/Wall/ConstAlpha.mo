within ThermalSeparation.Wall;
model ConstAlpha "constant alpha at column outside wall"
  extends ThermalSeparation.Wall.BaseWall;
  parameter SI.SurfaceCoefficientOfHeatTransfer alpha_const= 200
    "Heat transfer coefficient at the surface of the outer wall of the column";
equation
  for j in 1:n loop
    if geometry.hasInsulation then
      Qdot_out[j]=1/(1/alpha_const+geometry.s/geometry.lambda_ins) * 2  * Modelica.Constants.pi * geometry.H/n * (geometry.r2+geometry.r3)/2 * (T_wall[j]-T_amb);
    else
      Qdot_out[j]=alpha_const * 2  * Modelica.Constants.pi * geometry.H/n * (geometry.r2+geometry.r3)/2 * (T_wall[j]-T_amb);
    end if;
  end for;
  annotation (Icon(graphics),
                    Documentation(info=
          "<html>for information see <a href=\"Modelica://ThermalSeparation.Wall.BaseWall\">BaseWall</a></html>"));
end ConstAlpha;
