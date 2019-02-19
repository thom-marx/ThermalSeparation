within ThermalSeparation.Components.SourcesSinks;
model SinkGas
   extends Icons.Color.GasSink;
   replaceable package Medium=ThermalSeparation.Media.BaseMediumVapour
    "medium to be used"                                                                    annotation(choicesAllMatching);

  parameter Boolean use_p=true "set pressure to a constant value";
  parameter SI.Pressure p=1e5 "pressure in vapour sink"   annotation(Dialog(enable=use_p));
  ThermalSeparation.Interfaces.GasPortIn
                                gasPortIn(
                                 redeclare package Medium=Medium)       annotation (Placement(transformation(
          extent={{-96,0},{-76,20}},   rotation=0), iconTransformation(extent={{-116,
            -20},{-76,20}})));
 // Modelica.Blocks.Interfaces.RealInput Vdot_In
 //   annotation (Placement(transformation(extent={{-110,34},{-70,74}})));
equation
 if use_p then
  gasPortIn.p=p;
 end if;
  gasPortIn.x_outflow = fill(1,Medium.nSubstance);
  gasPortIn.h_outflow = 1;

  annotation (Documentation(info="<html>
<p>Vapour sink. By default the pressure in the vapour sink is fixed. If the pressure shall not be fixed, the parameter use_p has to be set to false.</p>
</html>"), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics));
end SinkGas;
