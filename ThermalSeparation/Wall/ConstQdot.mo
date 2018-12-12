within ThermalSeparation.Wall;
model ConstQdot "constant heat flow rate to ambiance"
  extends ThermalSeparation.Wall.BaseWall;
  parameter SI.HeatFlowRate Qdot_set = -500 "negative if heat is lost";
equation
  Qdot_out=fill(-Qdot_set/n, n);
  annotation (Icon(graphics),
Documentation(info=
          "<html>for information see <a href=\"Modelica://ThermalSeparation.Wall.BaseWall\">BaseWall</a></html>"));
end ConstQdot;
