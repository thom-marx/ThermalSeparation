within ThermalSeparation.Wall;
model SetQdot "real input for heat flow rate to ambiance"
  extends ThermalSeparation.Wall.BaseWall;
  parameter SI.HeatFlowRate Qdot_set = -500 "negative if heat is lost";
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-100,-10},{-60,30}})));
equation
  Qdot_out=fill(-u/n, n);
  annotation (Icon(graphics),
Documentation(info=
          "<html>for information see <a href=\"Modelica://ThermalSeparation.Wall.BaseWall\">BaseWall</a></html>"));
end SetQdot;
