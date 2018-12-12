within ThermalSeparation.Wall;
model Adiabatic "adiabatic: no heat loss to ambiance"
  extends ThermalSeparation.Wall.BaseWall;
equation
  Qdot_out=zeros(n);
  annotation (Icon(graphics),
Documentation(info=
          "<html>for information see <a href=\"Modelica://ThermalSeparation.Wall.BaseWall\">BaseWall</a></html>"));
end Adiabatic;
