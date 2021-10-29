within ThermalSeparation.Components.Compressors;
model FlueGasCompressor
  "simple model for flue gas compressor with fixed inlet pressure "
   extends Icons.Color.Compressors;
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
