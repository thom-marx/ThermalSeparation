within ThermalSeparation.Utilities;
block LimPID_Input
  "P, PI, PD, and PID controller with limited output, anti-windup compensation and setpoint weighting using inputs instead of connectors for u_s and u_m"
  import Modelica.Blocks.Types.SimpleController;
  import Modelica.Blocks.Types.Init;
  extends Modelica.Blocks.Icons.Block;

  Modelica.Blocks.Interfaces.RealOutput y "Connector of actuator output signal"
    annotation (Placement(transformation(extent={{100,-10},{120,10}},
        rotation=0)));
  input Real u_s;
  Real u_m;
  output Real controlError = u_s - u_m
    "Control error (set point - measurement)";

  parameter Modelica.Blocks.Types.SimpleController controllerType=
         Modelica.Blocks.Types.SimpleController.PID "Type of controller";
  parameter Real k(min=0) = 1 "Gain of controller";
  parameter Modelica.Units.SI.Time Ti(min=Modelica.Constants.small, start=0.5)=0.5
    "Time constant of Integrator block"
     annotation(Dialog(enable=controllerType==SimpleController.PI or
                              controllerType==SimpleController.PID));
  parameter Modelica.Units.SI.Time Td(min=0, start=0.1)=0.1
    "Time constant of Derivative block"
       annotation(Dialog(enable=controllerType==SimpleController.PD or
                                controllerType==SimpleController.PID));
  parameter Real yMax(start=1) "Upper limit of output";
  parameter Real yMin=-yMax "Lower limit of output";
  parameter Real wp(min=0) = 1 "Set-point weight for Proportional block (0..1)";
  parameter Real wd(min=0) = 0 "Set-point weight for Derivative block (0..1)"
       annotation(Dialog(enable=controllerType==SimpleController.PD or
                                controllerType==SimpleController.PID));
  parameter Real Ni(min=100*Modelica.Constants.eps) = 0.9
    "Ni*Ti is time constant of anti-windup compensation"
     annotation(Dialog(enable=controllerType==SimpleController.PI or
                              controllerType==SimpleController.PID));
  parameter Real Nd(min=100*Modelica.Constants.eps) = 10
    "The higher Nd, the more ideal the derivative block"
       annotation(Dialog(enable=controllerType==SimpleController.PD or
                                controllerType==SimpleController.PID));
  parameter Modelica.Blocks.Types.Init initType= Modelica.Blocks.Types.Init.InitialOutput
    "Type of initialization (1: no init, 2: steady state, 3: initial state, 4: initial output)"
                                     annotation(Evaluate=true,
      Dialog(group="Initialization"));
  parameter Boolean limitsAtInit = true
    "= false, if limits are ignored during initializiation"
    annotation(Evaluate=true, Dialog(group="Initialization",
                       enable=controllerType==SimpleController.PI or
                              controllerType==SimpleController.PID));
  parameter Real xi_start=0
    "Initial or guess value value for integrator output (= integrator state)"
    annotation (Dialog(group="Initialization",
                enable=controllerType==SimpleController.PI or
                       controllerType==SimpleController.PID));
  parameter Real xd_start=0
    "Initial or guess value for state of derivative block"
    annotation (Dialog(group="Initialization",
                         enable=controllerType==SimpleController.PD or
                                controllerType==SimpleController.PID));
  parameter Real y_start=0 "Initial value of output"
    annotation(Dialog(enable=initType == Init.InitialOutput, group=
          "Initialization"));
  constant SI.Time unitTime=1  annotation(HideResult=true);
  Modelica.Blocks.Math.Add addP(k1=wp, k2=-1)
    annotation (Placement(transformation(extent={{-80,40},{-60,60}}, rotation=
           0)));
  Modelica.Blocks.Math.Add addD(k1=wd, k2=-1) if with_D
    annotation (Placement(transformation(extent={{-80,-10},{-60,10}},
          rotation=0)));
  Modelica.Blocks.Math.Gain P(k=1)
                     annotation (Placement(transformation(extent={{-40,40},{
            -20,60}}, rotation=0)));
  Modelica.Blocks.Continuous.Integrator I(
    k=unitTime/Ti, y_start=xi_start,
    initType=if initType == Init.SteadyState then Init.SteadyState else if
        initType == Init.InitialState
         then Init.InitialState else Init.NoInit) if with_I
    annotation (Placement(transformation(extent={{-40,-60},{-20,-40}},
          rotation=0)));
  Modelica.Blocks.Continuous.Derivative D(
    k=Td/unitTime, T=max([Td/Nd, 1.e-14]), x_start=xd_start,
    initType=if initType == Init.SteadyState or initType == Init.InitialOutput
         then Init.SteadyState else if initType == Init.InitialState then
        Init.InitialState else Init.NoInit) if with_D
    annotation (Placement(transformation(extent={{-40,-10},{-20,10}},
          rotation=0)));
  Modelica.Blocks.Math.Gain gainPID(k=k)
                                annotation (Placement(transformation(extent={
            {30,-10},{50,10}}, rotation=0)));
  Modelica.Blocks.Math.Add3 addPID
                          annotation (Placement(transformation(
          extent={{0,-10},{20,10}}, rotation=0)));
  Modelica.Blocks.Math.Add3 addI(k2=-1) if with_I
                                         annotation (Placement(
        transformation(extent={{-80,-60},{-60,-40}}, rotation=0)));
  Modelica.Blocks.Math.Add addSat(k1=+1, k2=-1) if with_I
    annotation (Placement(transformation(
        origin={80,-50},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Modelica.Blocks.Math.Gain gainTrack(k=1/(k*Ni)) if with_I
    annotation (Placement(transformation(extent={{40,-80},{20,-60}}, rotation=
           0)));
  Modelica.Blocks.Nonlinear.Limiter limiter(
    uMax=yMax,
    uMin=yMin,
    homotopyType=Modelica.Blocks.Types.LimiterHomotopy.Linear)
    annotation (Placement(transformation(extent={{70,-10},{90,10}}, rotation=
            0)));
protected
  parameter Boolean with_I = controllerType==SimpleController.PI or
                             controllerType==SimpleController.PID annotation(Evaluate=true, HideResult=true);
  parameter Boolean with_D = controllerType==SimpleController.PD or
                             controllerType==SimpleController.PID annotation(Evaluate=true, HideResult=true);
public
  Modelica.Blocks.Sources.Constant Dzero(k=0) if not with_D
    annotation (Placement(transformation(extent={{-30,20},{-20,30}}, rotation=
           0)));
  Modelica.Blocks.Sources.Constant Izero(k=0) if not with_I
    annotation (Placement(transformation(extent={{10,-55},{0,-45}}, rotation=
            0)));
  Modelica.Blocks.Interfaces.RealOutput output_u_s
    "Connector of actuator output signal";
  Modelica.Blocks.Interfaces.RealOutput output_u_m
    "Connector of actuator output signal";
  Modelica.Blocks.Interfaces.RealInput input_u_m;
    Modelica.Blocks.Interfaces.RealInput input_u_s;
initial equation
 if initType==Init.InitialOutput then
    gainPID.y = y_start;
 end if;
equation

    output_u_s=u_s;
    output_u_m=u_m;
  assert(yMax >= yMin, "LimPID: Limits must be consistent. However, yMax (=" + String(yMax) +
                       ") < yMin (=" + String(yMin) + ")");
  if initType == Init.InitialOutput and (y_start < yMin or y_start > yMax) then
      Modelica.Utilities.Streams.error("LimPID: Start value y_start (=" + String(y_start) +
         ") is outside of the limits of yMin (=" + String(yMin) +") and yMax (=" + String(yMax) + ")");
  end if;
  assert(limitsAtInit or not limitsAtInit and y >= yMin and y <= yMax,
         "LimPID: During initialization the limits have been switched off.\n" +
         "After initialization, the output y (=" + String(y) +
         ") is outside of the limits of yMin (=" + String(yMin) +") and yMax (=" + String(yMax) + ")");

  connect(addP.y, P.u) annotation (Line(points={{-59,50},{-42,50}}, color={0,
          0,127}));
  connect(addD.y, D.u)
    annotation (Line(points={{-59,0},{-42,0}}, color={0,0,127}));
  connect(addI.y, I.u) annotation (Line(points={{-59,-50},{-42,-50}}, color={
          0,0,127}));
  connect(P.y, addPID.u1) annotation (Line(points={{-19,50},{-10,50},{-10,8},
          {-2,8}}, color={0,0,127}));
  connect(D.y, addPID.u2)
    annotation (Line(points={{-19,0},{-2,0}}, color={0,0,127}));
  connect(I.y, addPID.u3) annotation (Line(points={{-19,-50},{-10,-50},{-10,
          -8},{-2,-8}}, color={0,0,127}));
  connect(addPID.y, gainPID.u)
    annotation (Line(points={{21,0},{28,0}}, color={0,0,127}));
  connect(gainPID.y, addSat.u2) annotation (Line(points={{51,0},{60,0},{60,
          -20},{74,-20},{74,-38}}, color={0,0,127}));
  connect(gainPID.y, limiter.u)
    annotation (Line(points={{51,0},{68,0}}, color={0,0,127}));
  connect(limiter.y, addSat.u1) annotation (Line(points={{91,0},{94,0},{94,
          -20},{86,-20},{86,-38}}, color={0,0,127}));
  connect(limiter.y, y)
    annotation (Line(points={{91,0},{110,0}}, color={0,0,127}));
  connect(addSat.y, gainTrack.u) annotation (Line(points={{80,-61},{80,-70},{
          42,-70}}, color={0,0,127}));
  connect(gainTrack.y, addI.u3) annotation (Line(points={{19,-70},{-88,-70},{
          -88,-58},{-82,-58}}, color={0,0,127}));
  connect(Dzero.y, addPID.u2) annotation (Line(points={{-19.5,25},{-14,25},{
          -14,0},{-2,0}}, color={0,0,127}));
  connect(Izero.y, addPID.u3) annotation (Line(points={{-0.5,-50},{-10,-50},{
          -10,-8},{-2,-8}}, color={0,0,127}));
  connect(output_u_s, input_u_s)
                       annotation (Line(
      points={{-101.5,-11.5},{-101.5,-23.25},{-105,-23.25},{-105,-37}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(input_u_s, addD.u1)
                            annotation (Line(
      points={{-105,-37},{-105,-13.5},{-82,-13.5},{-82,6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(input_u_s, addI.u1)
                            annotation (Line(
      points={{-105,-37},{-96,-37},{-96,-42},{-82,-42}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(input_u_s, addP.u1)
                            annotation (Line(
      points={{-105,-37},{-105,11.5},{-82,11.5},{-82,56}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(output_u_m, input_u_m)
                        annotation (Line(
      points={{1.5,-82.5},{-25.75,-82.5},{-25.75,-89},{-48,-89}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(input_u_m, addI.u2)
                             annotation (Line(
      points={{-48,-89},{-48,-70},{-82,-70},{-82,-50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(input_u_m, addD.u2)
                             annotation (Line(
      points={{-48,-89},{-48,-47},{-82,-47},{-82,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(input_u_m, addP.u2)
                             annotation (Line(
      points={{-48,-89},{-48,-22.5},{-82,-22.5},{-82,44}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (defaultComponentName="PID",
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={1,1}), graphics={
        Line(points={{-80,78},{-80,-90}}, color={192,192,192}),
        Polygon(
          points={{-80,90},{-88,68},{-72,68},{-80,90}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,-80},{82,-80}}, color={192,192,192}),
        Polygon(
          points={{90,-80},{68,-72},{68,-88},{90,-80}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-80,-80},{-80,50},{-80,-20},{30,60},{80,60}}, color={0,
              0,127}),
        Text(
          extent={{-20,-20},{80,-60}},
          lineColor={192,192,192},
          textString="PID")}),
    Documentation(info="<html>
<p>This model is taken completely from <a href=\"Modelica://Modelica.Blocks.Continuous.LimPID\">Modelica.Blocks.Continuous.LimPID</a>. The only difference is, that the values for u_s (setpoint value) and u_m (measurement signal) are submitted to the model via input variables and not using connectors. This approach is more suitable for the use of this controller in <a href=\"Modelica://ThermalSeparation.Components.Columns.BaseColumn\">BaseColumn</a>, since it avoids unnecessary generation of warnings in Dymola.</p>
<p>The following rest of the documentation is also from <a href=\"Modelica://Modelica.Blocks.Continuous.LimPID\">Modelica.Blocks.Continuous.LimPID</a>: </p>
<p><br/>Via parameter <b>controllerType</b> either <b>P</b>, <b>PI</b>, <b>PD</b>, or <b>PID</b> can be selected. If, e.g., PI is selected, all components belonging to the D-part are removed from the block (via conditional declarations). The example model <a href=\"Modelica://Modelica.Blocks.Examples.PID_Controller\">Modelica.Blocks.Examples.PID_Controller</a> demonstrates the usage of this controller. Several practical aspects of PID controller design are incorporated according to chapter 3 of the book: </p>
<dl><dt>&Aring;str&ouml;m K.J., and H&auml;gglund T.:</dt>
<dd><b>PID Controllers: Theory, Design, and Tuning</b>. Instrument Society of America, 2nd edition, 1995. </dd>
</dl><p>Besides the additive <b>proportional, integral</b> and <b>derivative</b> part of this controller, the following features are present: </p>
<p><ul>
<li>The output of this controller is limited. If the controller is in its limits, anti-windup compensation is activated to drive the integrator state to zero. </li>
<li>The high-frequency gain of the derivative part is limited to avoid excessive amplification of measurement noise.</li>
<li>Setpoint weighting is present, which allows to weight the setpoint in the proportional and the derivative part independantly from the measurement. The controller will respond to load disturbances and measurement noise independantly of this setting (parameters wp, wd). However, setpoint changes will depend on this setting. For example, it is useful to set the setpoint weight wd for the derivative part to zero, if steps may occur in the setpoint signal. </li>
</ul></p>
<p>The parameters of the controller can be manually adjusted by performing simulations of the closed loop system (= controller + plant connected together) and using the following strategy: </p>
<p><ol>
<li>Set very large limits, e.g., yMax = Modelica.Constants.inf</li>
<li>Select a <b>P</b>-controller and manually enlarge parameter <b>k</b> (the total gain of the controller) until the closed-loop response cannot be improved any more.</li>
<li>Select a <b>PI</b>-controller and manually adjust parameters <b>k</b> and <b>Ti</b> (the time constant of the integrator). The first value of Ti can be selected, such that it is in the order of the time constant of the oscillations occuring with the P-controller. If, e.g., vibrations in the order of T=10 ms occur in the previous step, start with Ti=0.01 s.</li>
<li>If you want to make the reaction of the control loop faster (but probably less robust against disturbances and measurement noise) select a <b>PID</b>-Controller and manually adjust parameters <b>k</b>, <b>Ti</b>, <b>Td</b> (time constant of derivative block).</li>
<li>Set the limits yMax and yMin according to your specification.</li>
<li>Perform simulations such that the output of the PID controller goes in its limits. Tune <b>Ni</b> (Ni*Ti is the time constant of the anti-windup compensation) such that the input to the limiter block (= limiter.u) goes quickly enough back to its limits. If Ni is decreased, this happens faster. If Ni=infinity, the anti-windup compensation is switched off and the controller works bad. </li>
</ol></p>
<p><b>Initialization</b> </p>
<p>This block can be initialized in different ways controlled by parameter <b>initType</b>. The possible values of initType are defined in <a href=\"Modelica://Modelica.Blocks.Types.InitPID\">Modelica.Blocks.Types.InitPID</a>. This type is identical to <a href=\"Modelica://Modelica.Blocks.Types.Init\">Types.Init</a>, with the only exception that the additional option <b>DoNotUse_InitialIntegratorState</b> is added for backward compatibility reasons (= integrator is initialized with InitialState whereas differential part is initialized with NoInit which was the initialization in version 2.2 of the Modelica standard library). </p>
<p>Based on the setting of initType, the integrator (I) and derivative (D) blocks inside the PID controller are initialized according to the following table: </p>
<table cellspacing=\"0\" cellpadding=\"2\" border=\"1\"><tr>
<td valign=\"top\"><p><h4>initType</h4></p></td>
<td valign=\"top\"><p><h4>I.initType</h4></p></td>
<td valign=\"top\"><p><h4>D.initType</h4></p></td>
</tr>
<tr>
<td valign=\"top\"><p><h4>NoInit</h4></p></td>
<td valign=\"top\"><p>NoInit</p></td>
<td valign=\"top\"><p>NoInit</p></td>
</tr>
<tr>
<td valign=\"top\"><p><h4>SteadyState</h4></p></td>
<td valign=\"top\"><p>SteadyState</p></td>
<td valign=\"top\"><p>SteadyState</p></td>
</tr>
<tr>
<td valign=\"top\"><p><h4>InitialState</h4></p></td>
<td valign=\"top\"><p>InitialState</p></td>
<td valign=\"top\"><p>InitialState</p></td>
</tr>
<tr>
<td valign=\"top\"><p><h4>InitialOutput</h4></p><p>and initial equation: y = y_start</p></td>
<td valign=\"top\"><p>NoInit</p></td>
<td valign=\"top\"><p>SteadyState</p></td>
</tr>
<tr>
<td valign=\"top\"><p><h4>DoNotUse_InitialIntegratorState</h4></p></td>
<td valign=\"top\"><p>InitialState</p></td>
<td valign=\"top\"><p>NoInit</p></td>
</tr>
</table>
<p><br/><br/>In many cases, the most useful initial condition is <b>SteadyState</b> because initial transients are then no longer present. If initType = Init.SteadyState, then in some cases difficulties might occur. The reason is the equation of the integrator: </p>
<p><code><b>der</b></code><code>(y) = k*u;</code> </p>
<p>The steady state equation &QUOT;der(x)=0&QUOT; leads to the condition that the input u to the integrator is zero. If the input u is already (directly or indirectly) defined by another initial condition, then the initialization problem is <b>singular</b> (has none or infinitely many solutions). This situation occurs often for mechanical systems, where, e.g., u = desiredSpeed - measuredSpeed and since speed is both a state and a derivative, it is natural to initialize it with zero. As sketched this is, however, not possible. The solution is to not initialize u_m or the variable that is used to compute u_m by an algebraic equation. </p>
<p>If parameter <b>limitAtInit</b> = <b>false</b>, the limits at the output of this controller block are removed from the initialization problem which leads to a much simpler equation system. After initialization has been performed, it is checked via an assert whether the output is in the defined limits. For backward compatibility reasons <b>limitAtInit</b> = <b>true</b>. In most cases it is best to use <b>limitAtInit</b> = <b>false</b>. </p>
</html>"),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={1,1}), graphics),
    uses(Modelica(version="3.1")));
end LimPID_Input;
