within ThermalSeparation.Components.Compressors;
package OldCompressorModels

  model Compressor_GasPortInOut
    "Base class for components with two fluid portsTable look-up in two dimensions (matrix/file) "
       extends Icons.Icons.Compressors;
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

   parameter SI.Pressure paNom = 4.69
      "Nominal inlet pressure for which table is valid";

    Integer tableID;

    import SI = Modelica.SIunits;
  //  constant Real PI =  Modelica.Constants.pi;

    replaceable package Medium =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_N2
                                                             constrainedby ThermalSeparation.Media.BaseMediumVapour
                                                                                                        annotation(choicesAllMatching);

    parameter SI.Pressure p_start_in = Medium.reference_p
      "Start value of inlet pressure"
      annotation(Dialog(tab = "Initialization"));
    parameter SI.Pressure p_start_out = Medium.reference_p
      "Start value of outlet pressure"
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.Temperature T_start_in= Medium.reference_T
      "Start value of temperature"
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.Temperature T_start_out= Medium.reference_T
      "Start value of  outlet temperature"
      annotation(Dialog(tab = "Initialization", enable = use_T_start));

  //  SI.MoleFraction x_in[Medium.nS] = Medium.reference_X;

    parameter Modelica.SIunits.MassFlowRate m_start=15
      "Guess value for mass flow rate"
      annotation(Dialog(tab="Initialization"));
    parameter Medium.MassFraction X_start_in[Medium.nX] = Medium.reference_X
      "Start value of mass fractions m_i/m"
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter Medium.MassFraction X_start_out[Medium.nX] = Medium.reference_X
      "Start value of mass fractions m_i/m"
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    final parameter Integer nS=Medium.nSubstance "number of species";

   /* parameter Boolean preferredStates=true 
    "Try to select preferred medium states"                                      annotation(Dialog(tab="Advanced"));
// Initializatin parameters
//  constant Real k = 1.314;
*/
   parameter Real alpha = -20.;
   parameter Boolean use_alpha_in = false;

   Modelica.Blocks.Interfaces.RealInput alpha_in if use_alpha_in  annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-60,94}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-120,0})));
   Medium.SpecificEnthalpy deltah;

   Medium.ThermodynamicState stateIn(p(start=p_start_in), T(start=T_start_in))
      "state for medium inflowing through GasPortIn";
   Medium.ThermodynamicState stateOut(p(start=p_start_out), T(start=T_start_out))
      "state for medium inflowing through GasPortOut";
   Medium.AbsolutePressure deltaPressure;
    Interfaces.GasPortIn gasPortIn(redeclare package Medium=Medium)
      annotation (Placement(transformation(extent={{-20,40},{0,60}}),
          iconTransformation(extent={{-20,-114},{20,-74}})));
    Interfaces.GasPortOut gasPortOut(redeclare package Medium=Medium)
      annotation (Placement(transformation(extent={{-10,-100},{10,-80}}),
          iconTransformation(extent={{-24,74},{24,114}})));

  Real sum_x=sum(gasPortOut.x_outflow);

  equation
    connect(alpha_in, alpha_in_internal);
    if not use_alpha_in then
      alpha_in_internal = alpha;
    end if;

    //gasPortIn.c = gasPortOut.c;
    inStream(gasPortIn.x_outflow) = gasPortOut.x_outflow;
    gasPortIn.x_outflow=gasPortOut.x_outflow;
    //gasPortOut.p_medium= gasPortOut.p;

  //  for i in 1:nS loop
  //    gasPortIn.x[i] = gasPortIn.c[i] * Medium.molarMass(stateIn) /Medium.density(stateIn);
  //  end for;

    stateIn = Medium.setState_pTX(gasPortIn.p, gasPortIn.T, Medium.reference_X); //TODO inStream(gasPortIn.Xi_outflow)
    stateOut = Medium.setState_pTX(gasPortOut.p, gasPortOut.T, Medium.reference_X); //TODO inStream(gasPortOut.Xi_outflow)
   // mediumIn = Medium.setState_pTX(gasPortIn.p, gasPortIn.T, Medium.reference_X);
   // mediumOut = Medium.setState_pTX(gasPortOut.p, gasPortOut.T, Medium.reference_X);

    // mass balance
    gasPortIn.Vdot*Medium.density(stateIn) + gasPortOut.Vdot*Medium.density(stateOut) = 0.0;//Vdot*rho

    // energy balance
  //  deltah = Medium.isentropicEnthalpy(gasPortOut.p, stateIn)- inStream(gasPortIn.h_outflow);
    deltah = Medium.isentropicEnthalpy(gasPortOut.p, stateIn)- Medium.specificEnthalpy(stateIn);

   gasPortOut.h_outflow = inStream(gasPortIn.h_outflow) + deltah;
    Medium.specificEnthalpy(stateOut) = Medium.specificEnthalpy(stateIn) + deltah;
   gasPortIn.h_outflow = inStream(gasPortOut.h_outflow);

    if tableOnFile then
      assert(tableName<>"NoName", "tableOnFile = true and no table name given");
    end if;
    if not tableOnFile then
      assert(size(table,1) > 0 and size(table,2) > 0, "tableOnFile = false and parameter table is an empty matrix");
    end if;

    when initial() then
      tableID=tableInit(if tableOnFile then tableName else "NoName",
                        if tableOnFile then fileName else "NoName", table, smoothness);
    end when;

  //  gasPortOut.p - gasPortIn.p + paNom = tableIpo(tableID, alpha_in_internal, gasPortIn.m_flow)*1e5;
    deltaPressure = tableIpo(tableID, alpha_in_internal, gasPortIn.Vdot*Medium.density(stateIn))*1e5 - paNom;
    gasPortOut.p - gasPortIn.p = deltaPressure;

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
                     graphics), Diagram(coordinateSystem(preserveAspectRatio=
              false, extent={{-100,-100},{100,100}}), graphics));
  end Compressor_GasPortInOut;

  model Compressor_GasPortInOut_pconst
    "Base class for components with two fluid portsTable look-up in two dimensions (matrix/file); as CompressorGasPortInOut but considering isentropic efficiency and with a fixed pressure at the inlet "
       extends Icons.Icons.Compressors;
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

   parameter SI.Pressure paNom = 4.69
      "Nominal inlet pressure for which table is valid";

    Integer tableID;

    import SI = Modelica.SIunits;
  //  constant Real PI =  Modelica.Constants.pi;

    replaceable package Medium =
        ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_N2
                                                             constrainedby ThermalSeparation.Media.BaseMediumVapour
                                                                                                        annotation(choicesAllMatching);

    parameter SI.Pressure p_start_in = Medium.reference_p
      "Start value of inlet pressure"
      annotation(Dialog(tab = "Initialization"));
    parameter SI.Pressure p_start_out = Medium.reference_p
      "Start value of outlet pressure"
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.Temperature T_start_in= Medium.reference_T
      "Start value of temperature"
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.Temperature T_start_out= Medium.reference_T
      "Start value of  outlet temperature"
      annotation(Dialog(tab = "Initialization", enable = use_T_start));

  //  SI.MoleFraction x_in[Medium.nS] = Medium.reference_X;

    parameter Modelica.SIunits.MassFlowRate m_start=15
      "Guess value for mass flow rate"
      annotation(Dialog(tab="Initialization"));
    parameter Medium.MassFraction X_start_in[Medium.nX] = Medium.reference_X
      "Start value of mass fractions m_i/m"
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter Medium.MassFraction X_start_out[Medium.nX] = Medium.reference_X
      "Start value of mass fractions m_i/m"
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    final parameter Integer nS=Medium.nSubstance "number of species";

   /* parameter Boolean preferredStates=true 
    "Try to select preferred medium states"                                      annotation(Dialog(tab="Advanced"));
// Initializatin parameters
//  constant Real k = 1.314;
*/
   parameter Real alpha = -20.;
   parameter Boolean use_alpha_in = false;

   parameter Real eta_is = 1 "isentropic efficiency";

   Modelica.Blocks.Interfaces.RealInput alpha_in if use_alpha_in  annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-60,94}), iconTransformation(extent={{-140,-20},{-100,20}},
            rotation=0)));
   Medium.SpecificEnthalpy deltah;

   Medium.ThermodynamicState stateIn(p(start=p_start_in), T(start=T_start_in))
      "state for medium inflowing through GasPortIn";
   Medium.ThermodynamicState stateOut(p(start=p_start_out), T(start=T_start_out))
      "state for medium inflowing through GasPortOut";
   Medium.AbsolutePressure deltaPressure;
    Interfaces.GasPortIn gasPortIn(redeclare package Medium=Medium)
      annotation (Placement(transformation(extent={{0,-94},{20,-74}}),
          iconTransformation(extent={{-20,-114},{20,-74}})));
    Interfaces.GasPortOut gasPortOut(redeclare package Medium=Medium)
      annotation (Placement(transformation(extent={{0,94},{20,114}}),
          iconTransformation(extent={{-20,74},{20,114}})));

  Real sum_x=sum(gasPortOut.x);
  //Real xxx= alpha_in_internal;//if time < 100 then alpha_in_internal else min(20,max(-20,alpha_in_internal));
  parameter Real eta_el = 0.82; //angepasst, damit die Leistung für den 100%-Fall mit der Messung übereinstimmt
  SI.Power P_el = gasPortIn.Vdot*Medium.density(stateIn)* (Medium.specificEnthalpy(stateOut) - Medium.specificEnthalpy(stateIn))/eta_el;

  parameter SI.Pressure p_const_in = 0.96e5 "constant inlet pressure";
  equation
    connect(alpha_in, alpha_in_internal);
    if not use_alpha_in then
      alpha_in_internal = alpha;
    end if;

    gasPortOut.c =gasPortIn.c*Medium.density(stateIn)/Medium.density(stateOut);
    gasPortIn.x = gasPortOut.x;
    gasPortOut.p_medium= gasPortOut.p;

  //  for i in 1:nS loop
  //    gasPortIn.x[i] = gasPortIn.c[i] * Medium.molarMass(stateIn) /Medium.density(stateIn);
  //  end for;

    stateIn = Medium.setState_pTX(gasPortIn.p, gasPortIn.T, Medium.reference_X); //TODO inStream(gasPortIn.Xi_outflow)
    stateOut = Medium.setState_pTX(gasPortOut.p, gasPortOut.T, Medium.reference_X); //TODO inStream(gasPortOut.Xi_outflow)
   // mediumIn = Medium.setState_pTX(gasPortIn.p, gasPortIn.T, Medium.reference_X);
   // mediumOut = Medium.setState_pTX(gasPortOut.p, gasPortOut.T, Medium.reference_X);

    // mass balance
    gasPortIn.Vdot*Medium.density(stateIn) + gasPortOut.Vdot*Medium.density(stateOut) = 0.0;//Vdot*rho

    // energy balance
  //  deltah = Medium.isentropicEnthalpy(gasPortOut.p, stateIn)- inStream(gasPortIn.h_outflow);
    deltah = Medium.isentropicEnthalpy(gasPortOut.p, stateIn)- Medium.specificEnthalpy(stateIn);

  //  gasPortOut.h_outflow = inStream(gasPortIn.h_outflow) + deltah;
    Medium.specificEnthalpy(stateOut) = Medium.specificEnthalpy(stateIn) + deltah/eta_is;
  //  gasPortIn.h_outflow = inStream(gasPortOut.h_outflow);

    if tableOnFile then
      assert(tableName<>"NoName", "tableOnFile = true and no table name given");
    end if;
    if not tableOnFile then
      assert(size(table,1) > 0 and size(table,2) > 0, "tableOnFile = false and parameter table is an empty matrix");
    end if;

    when initial() then
      tableID=tableInit(if tableOnFile then tableName else "NoName",
                        if tableOnFile then fileName else "NoName", table, smoothness);
    end when;

  //  gasPortOut.p - gasPortIn.p + paNom = tableIpo(tableID, alpha_in_internal, gasPortIn.m_flow)*1e5;
    gasPortOut.p - gasPortIn.p = deltaPressure;

  deltaPressure = tableIpo(tableID, alpha_in_internal, gasPortIn.Vdot*Medium.density(stateIn))*1e5 - paNom;
    gasPortIn.p = p_const_in;

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
  end Compressor_GasPortInOut_pconst;

  model FlueGasCompressor
    "simple model for flue gas compressor with fixed inlet pressure "
     extends Icons.Icons.Compressors;
        outer ThermalSeparation.SystemTS systemTS;
      parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

    import SI = Modelica.SIunits;

     replaceable package Medium =
         ThermalSeparation.Media.IdealGasMixtures.H2O_O2_CO2_N2
                                                              constrainedby ThermalSeparation.Media.BaseMediumVapour
                                                                                                         annotation(choicesAllMatching);

    Medium.BaseProperties mediumIn(c=gasPortIn.c, T0=T_ref, p=gasPortIn.p, T=gasPortIn.T, x=gasPortIn.x,  x_star=gasPortIn.x);
    Medium.BaseProperties mediumOut(c=gasPortOut.c, T0=T_ref,p=gasPortOut.p, T=gasPortOut.T, x=gasPortOut.x,  x_star=gasPortOut.x);
    parameter SI.Pressure p_start_in = Medium.reference_p
      "Start value of inlet pressure"
      annotation(Dialog(tab = "Initialization"));
    parameter SI.Pressure p_start_out = Medium.reference_p
      "Start value of outlet pressure"
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.Temperature T_start_in= Medium.reference_T
      "Start value of temperature"
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.Temperature T_start_out= Medium.reference_T
      "Start value of  outlet temperature"
      annotation(Dialog(tab = "Initialization", enable = use_T_start));

  //  SI.MoleFraction x_in[Medium.nS] = Medium.reference_X;

    parameter Modelica.SIunits.MassFlowRate m_start=15
      "Guess value for mass flow rate"
      annotation(Dialog(tab="Initialization"));
    parameter Medium.MassFraction X_start_in[Medium.nX] = Medium.reference_X
      "Start value of mass fractions m_i/m"
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter Medium.MassFraction X_start_out[Medium.nX] = Medium.reference_X
      "Start value of mass fractions m_i/m"
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
   final parameter Integer nS=Medium.nSubstance "number of species";

   parameter Real eta_is = 1 "isentropic efficiency";

   Medium.ThermodynamicState stateIn(p(start=p_start_in), T(start=T_start_in))
      "state for medium inflowing through GasPortIn";
   Medium.ThermodynamicState stateOut(p(start=p_start_out), T(start=T_start_out))
      "state for medium inflowing through GasPortOut";
   Medium.AbsolutePressure deltaPressure;
    Interfaces.GasPortIn gasPortIn(redeclare package Medium=Medium)
      annotation (Placement(transformation(extent={{0,-94},{20,-74}}),
          iconTransformation(extent={{-20,-114},{20,-74}})));
    Interfaces.GasPortOut gasPortOut(redeclare package Medium =Medium)
      annotation (Placement(transformation(extent={{-14,22},{6,42}}),
          iconTransformation(extent={{-20,74},{20,114}})));

  Real sum_x=sum(gasPortOut.x);
  parameter SI.Pressure p_in = 1.01e5 "constant inlet pressure";

  parameter Real eta_el = 0.95;
  SI.Power P_el = gasPortIn.Vdot*Medium.density(stateIn)* (Medium.specificEnthalpy(stateOut) - Medium.specificEnthalpy(stateIn))/eta_el;
  SI.Temperature T_is "temperature after isentropic compression";
  equation

    gasPortOut.c =gasPortIn.c*mediumIn.d/mediumOut.d;
    gasPortIn.x = gasPortOut.x;
    gasPortOut.p_medium= gasPortOut.p;
    gasPortOut.T = gasPortIn.T + (T_is - gasPortIn.T)/eta_is;

    T_is = gasPortIn.T * (gasPortOut.p/gasPortIn.p)^(0.4/1.4);

    stateIn = Medium.setState_pTX(gasPortIn.p, gasPortIn.T, Medium.reference_X); //TODO inStream(gasPortIn.Xi_outflow)
    stateOut = Medium.setState_pTX(gasPortOut.p, gasPortOut.T, Medium.reference_X); //TODO inStream(gasPortOut.Xi_outflow)

    // mass balance
    gasPortIn.Vdot*mediumIn.d + gasPortOut.Vdot*mediumOut.d = 0.0;

   gasPortOut.p - gasPortIn.p = deltaPressure;

    gasPortIn.p = p_in;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                     graphics));
  end FlueGasCompressor;
end OldCompressorModels;
