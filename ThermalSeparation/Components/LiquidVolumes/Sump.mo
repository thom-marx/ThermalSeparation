within ThermalSeparation.Components.LiquidVolumes;
model Sump
 extends Icons.Color.Sump;
replaceable package MediumLiquid =
ThermalSeparation.Media.WaterBasedLiquid.CO2_H2O     constrainedby ThermalSeparation.Media.BaseMediumLiquid
                                                                                                      annotation(choicesAllMatching);

    replaceable package MediumVapour =
      ThermalSeparation.Media.IdealGasMixtures.H2O_CO2     constrainedby ThermalSeparation.Media.BaseMediumVapour
                                                                                                      annotation(choicesAllMatching);
  outer ThermalSeparation.SystemTS systemTS;
  parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

  MediumLiquid.BaseProperties mediumLiquid(T0=T_ref,p=p, T=T, x=x_l,h=h_l);
  MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref,p=portIn.p, T=T_l_in, x=x_l_in,h=h_l_in);
  MediumLiquid.BaseProperties mediumLiquidInRecirc(T0=T_ref,p=portInRecirc.p, T=T_l_in_recirc, x=x_l_in_recirc,h=h_l_in_recirc);

    parameter Boolean inertLiquid[nSL] = fill(false,nSL);
  parameter Integer nS=2
    "number of species which are equal in vapour and liquid phase";
  final parameter Integer nL=nSL-nS
    "number of additional substances which are only in liquid phase";

  final parameter Integer nSL = MediumLiquid.nSubstance;

/*** Medium properties ***/
  SI.Density rho_l_in = mediumLiquidIn.d;

  ThermalSeparation.Units.MolarEnthalpy h_l_in;//= mediumLiquidIn.h;
  SI.MolarMass MM_l_in = mediumLiquidIn.MM;

  SI.Density rho_l_in_recirc = mediumLiquidInRecirc.d;

  ThermalSeparation.Units.MolarEnthalpy h_l_in_recirc;//= mediumLiquidInRecirc.h;
  SI.MolarMass MM_l_in_recirc = mediumLiquidInRecirc.MM;

  SI.Density rho_l = mediumLiquid.d;
  SI.Concentration c_l[nSL];
  Real dummy1(stateSelect=StateSelect.prefer)=c_l[1];
  //Real dummy2(stateSelect=StateSelect.prefer)=c_l[2];
   ThermalSeparation.Units.MolarEnthalpy h_l;//= mediumLiquid.h;
   ThermalSeparation.Units.MolarEnthalpy u_l(stateSelect=StateSelect.prefer) = mediumLiquid.u;
   SI.MolarMass MM_l = mediumLiquid.MM;
  SI.MoleFraction x_l[nSL];

 SI.Temperature T;
  SI.Temperature T_l_in;
   SI.Temperature T_l_in_recirc;

 SI.VolumeFlowRate Vdot_l_out(start=1);

  SI.Pressure p(start=2e5);

  /*** geometry data ***/
    final parameter SI.Area A= Modelica.Constants.pi/4* d_volume^2;
  SI.Height level(stateSelect=StateSelect.default);

  parameter SI.Length height_in = -0.15;
  parameter SI.Length length_out = 0.15;
  parameter SI.Diameter d_volume = 0.025;
  parameter Real zeta=2;

  parameter SI.Height inletHeight = 0;
  parameter SI.Pressure p_ambient = 1e5;

/*** Initialization ***/
parameter SI.Temperature T_start=300 annotation(Dialog(tab="Initialization"));
parameter Real x_start[nSL] = {0.0226,1-0.0226} annotation(Dialog(tab="Initialization"));
parameter SI.Height level_start= 0.1 annotation(Dialog(tab="Initialization"));

SI.MoleFraction x_l_in[nSL] = inStream(portIn.x_outflow);
SI.MoleFraction x_l_in_recirc[nSL] = inStream(portInRecirc.x_outflow);

  ThermalSeparation.Interfaces.LiquidPortOut
                           portOut(redeclare package Medium=MediumLiquid)
    annotation (Placement(transformation(extent={{100,-100},{120,-80}}),
        iconTransformation(extent={{-20,-116},{20,-76}})));
  ThermalSeparation.Interfaces.LiquidPortIn
                          portInRecirc(redeclare package Medium=MediumLiquid)
    annotation (Placement(transformation(extent={{100,0},{120,20}}),
        iconTransformation(extent={{66,0},{106,40}})));

    ThermalSeparation.Interfaces.LiquidPortIn
                            portIn(redeclare package Medium=MediumLiquid)
    annotation (Placement(transformation(extent={{0,96},{20,116}}),
        iconTransformation(extent={{-20,76},{20,116}})));

  Modelica.Blocks.Interfaces.RealOutput y=level
    annotation (Placement(transformation(extent={{78,58},{98,78}}),
        iconTransformation(extent={{78,58},{98,78}})));

/*** monitoring ***/
SI.MassFlowRate mdot_out=Ndot_l_out*MM_l;
SI.Volume V_liq = A*level;
SI.MolarFlowRate Ndot_l_in= portIn.Ndot "total molar flow rate vapour";
SI.MolarFlowRate Ndot_l_in_recirc= portInRecirc.Ndot
    "total molar flow rate vapour";
SI.MolarFlowRate Ndot_l_out= Vdot_l_out*rho_l/MM_l
    "total molar flow rate vapour";

  ThermalSeparation.Interfaces.GasPortIn
                       gasPortIn(redeclare package Medium=MediumVapour)
    annotation (Placement(transformation(extent={{-160,-22},{-140,-2}}),
        iconTransformation(extent={{-104,14},{-64,58}})));
  ThermalSeparation.Interfaces.GasPortOut
                        gasPortOut(redeclare package Medium=MediumVapour)
    annotation (Placement(transformation(extent={{-160,0},{-140,20}}),
        iconTransformation(extent={{-104,58},{-64,98}})));
equation
      inStream(gasPortOut.h_outflow)=gasPortIn.h_outflow;
       inStream(gasPortOut.x_outflow)=gasPortIn.x_outflow;
         inStream(portOut.h_outflow)=portIn.h_outflow;
       inStream(portOut.x_outflow)=portIn.x_outflow;
            inStream(portOut.h_outflow)=portInRecirc.h_outflow;
       inStream(portOut.x_outflow)=portInRecirc.x_outflow;


    inStream(gasPortIn.h_outflow)=gasPortOut.h_outflow;
      inStream(gasPortIn.x_outflow)=gasPortOut.x_outflow;
        gasPortIn.p=gasPortOut.p;

          gasPortIn.Ndot+gasPortOut.Ndot=0;

  portOut.x_outflow = x_l;
  portOut.p = p;
  portOut.h_outflow = h_l;
  //portOut.Ndot = - max(0,Ndot_l_out);

  h_l_in = inStream(portIn.h_outflow);
    h_l_in_recirc = inStream(portInRecirc.h_outflow);

  /*** correlation between mole fraction x and concentration c ***/
   for i in 1:nSL loop
     x_l[i] = c_l[i] * MM_l /rho_l;
   end for;

     /***energy balance ***/
   A* der(sum(c_l[:])*u_l*level)  =  portIn.Ndot*h_l_in + portInRecirc.Ndot*h_l_in_recirc +portOut.Ndot*h_l;
     /*** mass balance ***/
    for i in 1:nSL loop
            A* der(c_l[i]*level) = portIn.Ndot*x_l_in[i] + portInRecirc.Ndot*x_l_in_recirc[i] +portOut.Ndot*x_l[i];
     end for;
    //   A* der(rho_l/MM_l*level) = Vdot_l_in*rho_l_in/MM_l_in + Vdot_l_in_recirc*rho_l_in_recirc/MM_l_in_recirc - Vdot_l_out*rho_l/MM_l;
    sum(x_l)=1;

  //p=level*9.8*rho_l+portIn.p- zeta* rho_l/2*((Vdot_l_out)/A)^2;
  Vdot_l_out =sqrt(max(0,-( p-portIn.p - level*9.8*rho_l)/(zeta*rho_l/2)))*A;

 //   p=level*9.8*rho_l+portIn.p- zeta* rho_l/2*((Vdot_l_in+Vdot_l_in_recirc)/A)^2;
portIn.p=portInRecirc.p;
portIn.p=portOut.p;

initial equation
   T=T_start;

  level=level_start;

  //x_l=x_start;
  for i in 1:nSL-1 loop
    x_l[i] = x_start[i];
    end for;

 // p=1.015e5;


























equation
  connect(gasPortIn, gasPortIn) annotation (Line(
      points={{-150,-12},{-150,-12},{-150,-12}},
      color={255,127,39},
      thickness=1,
      smooth=Smooth.None));
  annotation (Diagram(graphics), Icon(coordinateSystem(preserveAspectRatio=false,
                   extent={{-100,-100},{100,100}}),
                                      graphics),
    Documentation(info="<html>
<p>Model of a column sump. The sump is basically a liquid storage with variable liquid level. This is not a separation stage. The vapour phase is not modelled explicitly, the vapour inlet and outlet properties are identical. The vapour phase is only needed to supply the pressure which is exerted on the liquid.</p>
</html>"));
end Sump;
