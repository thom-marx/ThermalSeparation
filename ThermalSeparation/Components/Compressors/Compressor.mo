within ThermalSeparation.Components.Compressors;
model Compressor
  "Base class for components with two fluid portsTable look-up in two dimensions (matrix/file) "
   extends Icons.Color.Compressors;
  parameter Boolean tableOnFile=false
    "true, if table is defined on file or in function usertab"
    annotation(Dialog(group="table data definition"));
  parameter Real table[:, :]=fill(0.0,2,2)
    "table matrix (grid u1 = first column, grid u2 = first row; e.g. table=[0,0;0,1])"
       annotation(Dialog(group="table data definition", enable = not tableOnFile));
  parameter String tableName="NoName"
    "table name on file or in function usertab (see docu)"
       annotation(Dialog(group="table data definition", enable = tableOnFile));
  parameter String fileName="NoName" "file where matrix is stored"
       annotation(Dialog(group="table data definition", enable = tableOnFile,
                         __Dymola_loadSelector(filter="Text files (*.txt);;Matlab files (*.mat)",
                         caption="Open file in which table is present")));
  parameter Modelica.Blocks.Types.Smoothness smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments
    "smoothness of table interpolation"
  annotation(Dialog(group="table data interpretation"));

parameter Medium.AbsolutePressure paNom = 4.69
    "Nominal inlet pressure for which table is valid";
  Integer tableID;

  import SI = Modelica.SIunits;
  constant Real PI =  Modelica.Constants.pi;

  Medium.ThermodynamicState medium_a(p(start=p_start_in), T(start=T_start_in))
    "actual state at port_a";
  Medium.ThermodynamicState medium_b(p(start=p_start_out), T(start=T_start_out))
    "actual state at port_b";
// Medium
  replaceable package Medium =
      Modelica.Media.IdealGases.MixtureGases.SimpleNaturalGas
    constrainedby Modelica.Media.Interfaces.PartialMedium
           annotation (choicesAllMatching = true);

  parameter Boolean preferredStates=true
    "Try to select preferred medium states"                                      annotation(Dialog(tab="Advanced"));
// Initializatin parameters
  parameter Modelica.SIunits.MassFlowRate m_start=15
    "Guess value for mass flow rate"
    annotation(Dialog(tab="Initialization"));
  parameter Medium.AbsolutePressure p_start_in = Medium.reference_p
    "Start value of inlet pressure"
    annotation(Dialog(tab = "Initialization"));
  parameter Medium.AbsolutePressure p_start_out = Medium.reference_p
    "Start value of outlet pressure"
    annotation(Dialog(tab = "Initialization"));
  parameter Boolean use_T_start = false
    "Use T_start if true, otherwise h_start"
    annotation(Dialog(tab = "Initialization"), Evaluate=true);
                                       //changed
  parameter Medium.SpecificEnthalpy h_start_in=
    if use_T_start then Medium.specificEnthalpy_pTX(p_start_in, T_start_in,X_start_in) else Medium.h_default
    "Start value of specific enthalpy"
    annotation(Dialog(tab = "Initialization", enable = not use_T_start));
  parameter Medium.SpecificEnthalpy h_start_out=
    if use_T_start then Medium.specificEnthalpy_pTX(p_start_out, T_start_out,X_start_out) else Medium.h_default
    "Start value of specific outlet enthalpy"
    annotation(Dialog(tab = "Initialization", enable = not use_T_start));
  parameter Medium.Temperature T_start_in=
    if use_T_start then Medium.reference_T else Medium.temperature_phX(p_start_in,h_start_in,X_start_in)
    "Start value of temperature"
    annotation(Dialog(tab = "Initialization", enable = use_T_start));
  parameter Medium.Temperature T_start_out=
    if use_T_start then Medium.reference_T else Medium.temperature_phX(p_start_out,h_start_out,X_start_out)
    "Start value of  outlet temperature"
    annotation(Dialog(tab = "Initialization", enable = use_T_start));
  parameter Medium.MassFraction X_start_in[Medium.nX] = Medium.reference_X
    "Start value of mass fractions m_i/m"
    annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
  parameter Medium.MassFraction X_start_out[Medium.nX] = Medium.reference_X
    "Start value of mass fractions m_i/m"
    annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));

  parameter Real alpha = -20.;
//  constant Real k = 1.314;
  Medium.SpecificEnthalpy deltah;
  parameter Boolean use_alpha_in = false;
  Modelica.Blocks.Interfaces.RealInput alpha_in if use_alpha_in  annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-60,94}), iconTransformation(extent={{-140,-20},{-100,20}},
          rotation=0)));

  Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium, m_flow(start=m_start), h_outflow(start=h_start_in), p(start=p_start_in),Xi_outflow(start=X_start_in[1:Medium.nXi]))
    "Inlet port" annotation (Placement(transformation(extent={{-16,80},{14,110}},
                  rotation=0), iconTransformation(extent={{-16,80},{14,110}})));
  Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium, m_flow(start=-m_start), h_outflow(start=h_start_out), p(start=p_start_out),Xi_outflow(start=X_start_out[1:Medium.nXi]))
    "Outlet port" annotation (Placement(transformation(extent={{14,-106},{-16,
            -76}},
          rotation=0), iconTransformation(extent={{14,-106},{-16,-76}})));

  Medium.ThermodynamicState state_from_a(p(start=p_start_in), T(start=T_start_in))
    "state for medium inflowing through port_a";
  Medium.ThermodynamicState state_from_b(p(start=p_start_out), T(start=T_start_out))
    "state for medium inflowing through port_b";

equation
      connect(alpha_in, alpha_in_internal);
  if not use_alpha_in then
    alpha_in_internal = alpha;
  end if;

// medium states
  state_from_a = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow));
  state_from_b = Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow));
  medium_a = Medium.setState_phX(port_a.p, actualStream(port_a.h_outflow), actualStream(port_a.Xi_outflow));
  medium_b = Medium.setState_phX(port_b.p, actualStream(port_b.h_outflow), actualStream(port_b.Xi_outflow));

  // no substance storage
  port_a.Xi_outflow = inStream(port_b.Xi_outflow);
  port_b.Xi_outflow = inStream(port_a.Xi_outflow);

  port_a.C_outflow = inStream(port_b.C_outflow);
  port_b.C_outflow = inStream(port_a.C_outflow);

  // mass balance
  port_a.m_flow + port_b.m_flow = 0.0;//Vdot*rho

  // energy balance
  // delta h
  // medium_b.T = (medium_a.T) *(port_b.p/port_a.p)^((k-1)/k);
  deltah = Medium.isentropicEnthalpy(port_b.p, state_from_a)- inStream(port_a.h_outflow);
  port_b.h_outflow = inStream(port_a.h_outflow) + deltah;
  port_a.h_outflow = inStream(port_b.h_outflow);

  if tableOnFile then
    assert(tableName<>"NoName", "tableOnFile = true and no table name given");
  end if;
  if not tableOnFile then
    assert(size(table,1) > 0 and size(table,2) > 0, "tableOnFile = false and parameter table is an empty matrix");
  end if;

  port_b.p - port_a.p + paNom = tableIpo(tableID, alpha_in_internal, port_a.m_flow)*1e5;
  when initial() then
    tableID=tableInit(if tableOnFile then tableName else "NoName",
                      if tableOnFile then fileName else "NoName", table, smoothness);
  end when;

protected
  Modelica.Blocks.Interfaces.RealInput alpha_in_internal
    "Needed to connect to conditional connector";
  function tableInit
    "Initialize 2-dim. table defined by matrix (for details see: Modelica/C-Sources/ModelicaTables.h)"

    input String tableName;
    input String fileName;
    input Real table[ :, :];
    input Modelica.Blocks.Types.Smoothness smoothness;
    output Integer tableID;
  external "C" tableID = ModelicaTables_CombiTable2D_init(
                 tableName, fileName, table, size(table, 1), size(table, 2),
                 smoothness);
    annotation(Library="ModelicaExternalC");
  end tableInit;

  function tableIpo
    "Interpolate 2-dim. table defined by matrix (for details see: Modelica/C-Sources/ModelicaTables.h)"
    input Integer tableID;
    input Real u1;
    input Real u2;
    output Real value;
  external "C" value=ModelicaTables_CombiTable2D_interpolate(tableID, u1, u2);
    annotation(Library="ModelicaExternalC");
  end tableIpo;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                   graphics));
end Compressor;
