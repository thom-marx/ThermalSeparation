within ThermalSeparation.Components.Compressors;
model CompressorSimple "isentropic compressor"
   extends Icons.Color.Compressors;
replaceable package Medium =
    ThermalSeparation.Media.BaseMediumVapour                        constrainedby ThermalSeparation.Media.BaseMediumVapour annotation(choicesAllMatching=true);

Medium.BaseProperties mediumIn(p=gasPortIn.p,x=inStream(gasPortIn.x_outflow),x_star=inStream(gasPortIn.x_outflow),c=c,h=inStream(gasPortIn.h_outflow));
Medium.BaseProperties mediumOut(p=gasPortOut.p,x=gasPortOut.x_outflow,x_star=gasPortOut.x_outflow,c=c);

parameter Real isExp=1.3 "isentropic exponent";
parameter Modelica.Units.SI.Power P_drive_const=10000 "power input of fan motor" annotation(Dialog(enable=not useP));
parameter Boolean useP=false "use P_drive from real input";

// parameter Real charLine[:,:]=fill(
//       0.0,
//       0,
//       2) "characteristic line for dp/Vdot in (bar/m^3/s)";
// parameter Boolean startConstant=false "start with constant dp" annotation(Dialog(tab="Initialization"));
// parameter Boolean useHomotopy=false "use homotopy method for init of pressure" annotation(Dialog(tab="Initialization"));
// parameter Modelica.Units.SI.Pressure dp_nom=0.05*1e5
//     "nominal pressure difference"                                                 annotation(Dialog(enable=startConstant or useHomotopy,tab="Initialization"));
// parameter Real omega_k=0.05
//     "large value if change between constant and variable shall be steep" annotation(Dialog(enable=startConstant,tab="Initialization"));
// parameter Real omega_time = 50 "Wendepunkt der tanh-Funktion" annotation(Dialog(enable=startConstant,tab="Initialization"));
//    Real omega = 0.5*tanh(omega_k*(time-omega_time))+0.5;

Modelica.Units.SI.Power P_hyd;
Modelica.Units.SI.Power P_drive;
Modelica.Units.SI.Pressure deltap;
Modelica.Units.SI.Concentration c[Medium.nSubstance];
Modelica.Units.SI.Temperature T_in=mediumIn.T;
Modelica.Units.SI.Temperature T_out=mediumOut.T;

Modelica.Units.SI.VolumeFlowRate Vdot_in=gasPortIn.Ndot*mediumIn.MM/mediumIn.d;
Modelica.Units.SI.VolumeFlowRate Vdot_out=-gasPortOut.Ndot*mediumOut.MM/mediumOut.d;

  Interfaces.GasPortIn gasPortIn(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{0,-94},{20,-74}}),
        iconTransformation(extent={{-20,-114},{20,-74}})));
  Interfaces.GasPortOut gasPortOut(redeclare package Medium=Medium)
    annotation (Placement(transformation(extent={{-26,44},{-6,64}}),
        iconTransformation(extent={{-22,74},{22,114}})));
//   Modelica.Blocks.Tables.CombiTable1Ds combiTable1Ds(table=charLine, u=Vdot_in)
//     annotation (Placement(transformation(extent={{-20,-10},{0,10}})));

  Modelica.Blocks.Interfaces.RealInput u if useP annotation (Placement(transformation(
          extent={{-140,-20},{-100,20}}),iconTransformation(extent={{-140,-20},
            {-100,20}})));

protected
Modelica.Blocks.Interfaces.RealInput P_int;

equation
  if not useP then
    P_int = P_drive_const;
  end if;

  if useP then
    P_drive=P_int;
  else
    P_drive=P_drive_const;
  end if;

  // pressure drop/flow characteristics
  P_hyd=P_drive;
  Vdot_in=P_hyd/deltap;

//   // pressure drop using characteristic line
//   if startConstant then
//    if useHomotopy then
//      deltap = homotopy(actual=combiTable1Ds.y[1]*1e5*omega+dp_nom*(1-omega), simplified=dp_nom);
//    else
//      deltap = combiTable1Ds.y[1]*1e5*omega+dp_nom*(1-omega);
//    end if;
//   else
//    if useHomotopy then
//      deltap = homotopy(actual=combiTable1Ds.y[1]*1e5, simplified=dp_nom);
//    else
//      deltap = combiTable1Ds.y[1]*1e5;
//    end if;
//   end if;

  // converting x to c for medium model
  c=inStream(gasPortIn.x_outflow)*mediumIn.d/mediumIn.MM;

  // energy balance (ideal gas assumption)
  T_out/T_in=(gasPortOut.p/gasPortIn.p)^((isExp-1)/isExp);

  // mole balance
  0 = gasPortIn.Ndot+gasPortOut.Ndot;

  // pressure
  gasPortOut.p=gasPortIn.p+deltap;

  // ports
  gasPortIn.h_outflow=inStream(gasPortOut.h_outflow);
  gasPortIn.x_outflow=inStream(gasPortOut.x_outflow);
  gasPortOut.h_outflow=mediumOut.h;
  gasPortOut.x_outflow=inStream(gasPortIn.x_outflow);

connect(P_int,u);
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                         graphics));
end CompressorSimple;
