within ThermalSeparation.Utilities;
model Adapt_k "model to adapt k"
parameter Integer samplePeriod=30;
parameter Real k_start;
parameter Integer startTime=10;
Real error;

  Modelica.Blocks.Math.Add add_vap
    annotation (Placement(transformation(extent={{6,-6},{26,14}})));
  Modelica.Blocks.Tables.CombiTable1D table(       u={error}, table=[0,-1; 1e-7,
        -0.5; 1e-6,0; 1e-5,0.5; 1e-4,0.7; 1e-3,0.9; 1e-2,1; 1e-1,1.5])
    annotation (Placement(transformation(extent={{-72,40},{-52,60}})));
  Modelica.Blocks.Math.Product product_vap
    annotation (Placement(transformation(extent={{-32,0},{-12,20}})));
  Modelica.Blocks.Discrete.ZeroOrderHold zeroOrderHold(
    samplePeriod=samplePeriod,
    startTime=startTime,
    ySample(start=k_start))
                     annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={30,-32})));
  Modelica.Blocks.Discrete.ZeroOrderHold zeroOrderHold3(
    ySample(start=0),
    samplePeriod=samplePeriod,
    startTime=startTime)
                     annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-48,-30})));
  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{90,-6},{110,14}})));
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-120,30},{-80,70}})));
equation
  connect(add_vap.y, y) annotation (Line(
      points={{27,4},{100,4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(y, product_vap.u1) annotation (Line(
      points={{100,4},{100,58},{-40,58},{-40,16},{-34,16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(zeroOrderHold.u, y)     annotation (Line(
      points={{42,-32},{100,-32},{100,4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(zeroOrderHold.y,add_vap. u2) annotation (Line(
      points={{19,-32},{8,-32},{8,-26},{-10,-26},{-10,-2},{4,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(table.y[1], zeroOrderHold3.u)     annotation (Line(
      points={{-51,50},{-48,50},{-48,-4},{-74,-4},{-74,-16},{-60,-16},{-60,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(zeroOrderHold3.y,product_vap. u2) annotation (Line(
      points={{-37,-30},{-36,-30},{-36,4},{-34,4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(product_vap.y,add_vap. u1) annotation (Line(
      points={{-11,10},{4,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(u, table.u[1]) annotation (Line(
      points={{-100,50},{-74,50}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics));
end Adapt_k;
