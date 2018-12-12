within ThermalSeparation.Utilities;
block PI_Input
  "Proportional-Integral controller using input instead of connector"
  import Modelica.Blocks.Types.Init;
  parameter Real k=1 "Gain";
  parameter SI.Time T(start=1, min=Modelica.Constants.small)
    "Time Constant (T>0 required)";
  parameter Modelica.Blocks.Types.Init initType=Modelica.Blocks.Types.Init.NoInit
    "Type of initialization (1: no init, 2: steady state, 3: initial state, 4: initial output)"
                                                                            annotation(Evaluate=true,
      Dialog(group="Initialization"));
  parameter Real x_start=0 "Initial or guess value of state"
    annotation (Dialog(group="Initialization"));
  parameter Real y_start=0 "Initial value of output"
    annotation(Dialog(enable=initType == Init.SteadyState or initType == Init.InitialOutput, group=
          "Initialization"));

  extends Modelica.Blocks.Interfaces.BlockIcon;

  input Real u "Connector of Real input signal";
  Modelica.Blocks.Interfaces.RealOutput y "Connector of Real output signal"
    annotation (Placement(transformation(extent={{100,-10},{120,10}},
        rotation=0)));
  output Real x(start=x_start) "State of block";

initial equation
  if initType == Init.SteadyState then
    der(x) = 0;
  elseif initType == Init.InitialState then
    x = x_start;
  elseif initType == Init.InitialOutput then
    y = y_start;
  end if;
equation
  der(x) = u/T;
  y = k*(x + u);
  annotation (defaultComponentName="PI",
    Documentation(info="
<HTML>
<p>This model is taken completely from <a href=\"Modelica://Modelica.Blocks.Continuous.PI\">Modelica.Blocks.Continuous.PI</a>. The only difference is, that the values for u_s (setpoint value) and u_m (measurement signal) are submitted to the model via input variables and not using connectors. This approach is more suitable for the use of this controller in <a href=\"Modelica://ThermalSeparation.Components.Columns.BaseColumn\">BaseColumn</a>, since it avoids unnecessary generation of warnings in Dymola.</p>
<p>The following rest of the documentation is also from <a href=\"Modelica://Modelica.Blocks.Continuous.PI\">Modelica.Blocks.Continuous.PI</a>: </p>
<p>
This blocks defines the transfer function between the input u and
the output y (element-wise) as <i>PI</i> system:
</p>
<pre>
                 1
   y = k * (1 + ---) * u
                T*s
           T*s + 1
     = k * ------- * u
             T*s
</pre>
<p>
If you would like to be able to change easily between different
transfer functions (FirstOrder, SecondOrder, ... ) by changing
parameters, use the general model class <b>TransferFunction</b>
instead and model a PI SISO system with parameters<br>
b = {k*T, k}, a = {T, 0}.
</p>
<pre>
Example:

   parameter: k = 0.3,  T = 0.4

   results in:
               0.4 s + 1
      y = 0.3 ----------- * u
                 0.4 s
</pre>

<p>
It might be difficult to initialize the PI component in steady state
due to the integrator part.
This is discussed in the description of package
<a href=\"Modelica://Modelica.Blocks.Continuous#info\">Continuous</a>.
</p>

</HTML>
"), Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Line(points={{-80,78},{-80,-90}}, color={192,192,192}),
        Polygon(
          points={{-80,90},{-88,68},{-72,68},{-80,88},{-80,90}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,-80},{82,-80}}, color={192,192,192}),
        Polygon(
          points={{90,-80},{68,-72},{68,-88},{90,-80}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-80,-80},{-80,-20},{60,80}}, color={0,0,127}),
        Text(
          extent={{0,6},{60,-56}},
          lineColor={192,192,192},
          textString="PI"),
        Text(
          extent={{-150,-150},{150,-110}},
          lineColor={0,0,0},
          textString="T=%T")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Rectangle(extent={{-60,60},{60,-60}}, lineColor={0,0,255}),
        Text(
          extent={{-68,24},{-24,-18}},
          lineColor={0,0,0},
          textString="k"),
        Text(
          extent={{-32,48},{60,0}},
          lineColor={0,0,0},
          textString="T s + 1"),
        Text(
          extent={{-30,-8},{52,-40}},
          lineColor={0,0,0},
          textString="T s"),
        Line(points={{-24,0},{54,0}}, color={0,0,0}),
        Line(points={{-100,0},{-60,0}}, color={0,0,255}),
        Line(points={{62,0},{100,0}}, color={0,0,255})}));
end PI_Input;
