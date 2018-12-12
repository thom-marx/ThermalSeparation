within ThermalSeparation.Components.ColumnSump;
model Sump_Steam_1_2

    replaceable package MediumVapour = 
      ThermalSeparation.Media.C2H5OH_Water_Vap  constrainedby
    ThermalSeparation.Media.BaseMediumVapour                                             annotation(choicesAllMatching);

    MediumVapour.BaseProperties mediumVapour(p=p_system,T=T_system,x = x_v,x_star = x_v,c = c_v);

  replaceable package MediumLiquid = 
  ThermalSeparation.Media.C2H5OH_Water_Liq  constrainedby
    ThermalSeparation.Media.BaseMediumLiquid                                             annotation(choicesAllMatching);

    MediumLiquid.BaseProperties mediumLiquid(p=p_system,T=T_system,x=x_l);
      MediumLiquid.BaseProperties mediumLiquidIn(p=p_system,T=T_in,x=x_l_in);

parameter Integer nS=2;
parameter Integer nL=0;
parameter Integer nV=0;
  final parameter Integer nSL = nS + nL;
  final parameter Integer nSV = nS + nV;

parameter Boolean inertVapour[nSV] = fill(false,nSV)
    "true for each component which is inert";
parameter Boolean inertLiquid[nSL] = fill(false,nSL)
    "true for each component which is inert";

    // vorläufig konstant
//parameter Real gamma=1;
//parameter Real phi=1;
//parameter Real phi_sat=1;

parameter SI.HeatFlowRate Qdot_out=0;
SI.HeatFlowRate Qdot_in;

//inflow
     SI.VolumeFlowRate Vdot_in_l; // konstant
     //input SI.Pressure p_in;//=2e5;
     SI.MoleFraction x_l_in[nSL]=liquidIn.x;//={0.1,0.9};
     SI.Temperature T_in=liquidIn.T;//=353.15;
     SI.Concentration c_l_in[nSL]=liquidIn.c;
     ThermalSeparation.Units.MolarEnthalpy h_l_in;
     SI.Density rho_l_in;
     SI.MolarMass MM_l_in;
     SI.Pressure p_sat[nS];//={1e5,0.5e5};

     constant SI.Acceleration g=Modelica.Constants.g_n;
     Real nmol[nS]={1,1};

   //  parameter SI.HeatFlowRate Qdot_in=800;

/***BASE PROPERTIES***/
//Liquid
     SI.VolumeFlowRate Vdot_out_l;
     SI.Density rho_l;
     SI.Concentration c_l[nSL];
     SI.MoleFraction x_l[nSL](start={0.1,0.9});
     ThermalSeparation.Units.MolarEnthalpy h_l;
     ThermalSeparation.Units.MolarEnthalpy u_l;
     SI.MolarMass MM_l;
     SI.Temperature T_l(stateSelect=StateSelect.default);
     SI.Pressure p_out_l_port;
     SI.AmountOfSubstance HU_l;

//Vapour
     SI.Density rho_v;
     SI.Concentration c_v[nSV];
     SI.Concentration dummy(stateSelect=StateSelect.default)=c_v[2];
     SI.MoleFraction x_v[nSV](start={0.9,0.1});
     ThermalSeparation.Units.MolarEnthalpy h_v;
     ThermalSeparation.Units.MolarEnthalpy u_v;
     SI.MolarMass MM_v;
     SI.Temperature T_v(start=273.15);
     SI.VolumeFlowRate Vdot_out_v(start=0.25);
     SI.Pressure p_out_v_port;
     //parameter Real HU_v_start=0.1;
     SI.AmountOfSubstance HU_v;//(start=HU_v_start);

//    SI.VolumeFlowRate Vdot_v( start=1);
//    SI.VolumeFlowRate Vdot_l;

//Geometry
parameter SI.Area A=1;
parameter SI.Height h=2;
parameter SI.Length l_w=0.8; // Wehrlaenge
SI.Volume Vges;
SI.Height height_l(stateSelect=StateSelect.default);
SI.Height h_w; //weir height

SI.Pressure p_system(start=1.08e5);
SI.Temperature T_system;
 SI.Volume V_l;
 SI.Volume V_v;

Real N_H2O = c_l[2]*V_l + c_v[2]*V_v;
Real N_Eth =  c_l[1]*V_l + c_v[1]*V_v;

      Interfaces.LiquidPortIn liquidIn(nS=nSL,Vdot=Vdot_in_l) 
        annotation (Placement(transformation(extent={{8,60},{28,80}}),
        iconTransformation(extent={{8,60},{28,80}})));
     Interfaces.LiquidPortOut liquidOut( nS=nSL,p=p_out_l_port,x=x_l,c=c_l,T=T_l,Vdot=-Vdot_out_l) 
       annotation (Placement(transformation(extent={{0,-68},{20,-48}}),
        iconTransformation(extent={{-22,-66},{-2,-46}})));

     Interfaces.GasPortOut vapourOut( nS=nSV,x=x_v,c=c_v,T=T_v,Vdot=-Vdot_out_v,p=p_system, p_medium=p_system) 
       annotation (Placement(transformation(extent={{30,60},{50,80}}),
        iconTransformation(extent={{-52,60},{-32,80}})));

/* heat transfer */
parameter SI.Pressure p_steam=2e5;
parameter SI.VolumeFlowRate vdot_steam_in=1;

replaceable model HeatTransferModel = 
ThermalSeparation.Components.ColumnSump.HeatTransferModel.dT_m_2 constrainedby
    ThermalSeparation.Components.ColumnSump.HeatTransferModel.BaseHeatTransferModel
                                                                             annotation(choicesAllMatching=true);

HeatTransferModel heatTransferModel(p_s=p_steam,vdot_v_in=vdot_steam_in,T_f_in=T_in,T_f_out=T_system,p_sump=p_system,x_l=x_l_in,mdot_l_in=Vdot_in_l*rho_l_in);

/* control of liquid level */
  Modelica.Blocks.Continuous.LimPID PID(
    yMax=1,
    yMin=0.01,
    y_start=1,
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    k=10,
    Ti=0.8)                                               annotation (Placement(transformation(extent={{-26,20},{-16,30}})));
  Modelica.Blocks.Sources.RealExpression Setpoint(y=1);

/* control parameters */
Real z=heatTransferModel.vdot_v_out/vdot_steam_in;

initial equation
  T_l=95+273.15;
  height_l=1;//1.901;
  //Vdot_out_v=0.25;

equation

//Medienübergabe

   MM_v = mediumVapour.MM;//max(1e-6,sum(x_v_in.*MediumVapour.MMX));
   rho_v = mediumVapour.d;
   h_v = mediumVapour.h;
   u_v = mediumVapour.u;

   h_l_in =mediumLiquidIn.h;
   MM_l_in=mediumLiquidIn.MM;
   rho_l_in=mediumLiquidIn.d;
   //mediumLiquidIn.n_mol=nmol;

  //mediumLiquid.n_mol=nmol;
  h_l =mediumLiquid.h;
  u_l =mediumLiquid.u;
  rho_l = mediumLiquid.d;
  MM_l = mediumLiquid.MM;
  for j in 1:nS loop
  p_sat[j]=mediumLiquid.p_sat[j];
  end for;

  Setpoint.y=PID.u_s;
  PID.u_m=height_l;
  PID.y= h_w;

/*** correlation between mole fraction x and concentration c ***/
  for i in 1:nS loop
    x_l[i] = c_l[i] * MM_l /rho_l;
    x_v[i] = c_v[i] * MM_v /rho_v;
    //x_l_in[i] = c_l_in[i] * MM_l_in /rho_l_in;
  end for;
//   for i in 1:nL loop
//     x_l[i+nS] = c_l[i+nS] * MM_l /rho_l;
//     x_l_in[i+nS] = c_l_in[i+nS] * MM_l_in /rho_l_in;
//   end for;
//   for i in 1:nV loop
//     x_v[i+nS] = c_v[i+nS] * MM_v /rho_v;
//   end for;

//Geometry
  Vges=A*h;
  Vdot_out_l*rho_l/MM_l=if height_l > h_w then rho_l/MM_l*1.848*l_w*(abs(height_l-h_w))^1.5 else 0;
//   A*der(height_l)=Vdot_in_l-Vdot_out_l;
  A*height_l=HU_l*MM_l/rho_l;

  //A*der(height_l)+A*der(h-height_l)=Vdot_in_l-Vdot_out_l-Vdot_out_v;

  Vges=HU_l*MM_l/rho_l+HU_v*MM_v/rho_v;

//MoleBalance
    for j in 1:nS loop
      der(HU_v*x_v[j]+HU_l*x_l[j])=Vdot_in_l*rho_l_in/MM_l_in*x_l_in[j]-Vdot_out_l*rho_l/MM_l*x_l[j]-Vdot_out_v*rho_v/MM_v*x_v[j];//+v[j]*r /*Reaktion
    end for;

//      for i in 1:nSV loop
//      if inertVapour[i] then
//        der(HU_v*x_v[i]) =  - Vdot_out_v * x_v[i];
//      end if;
//      end for;
//      for i in 1:nSL loop
//      if inertLiquid[i] then
//        der(HU_l*x_l[i]) = Vdot_in_l*x_l_in[i]  - Vdot_out_l * x_l[i];
//      end if;
//      end for;

//Phase equilibrium
    for j in 1:nS loop
    x_l[j]*p_sat[j]=x_v[j]*p_system;
    end for;
    sum(x_l)=1;
    sum(x_v)=1;

//Enthalpy Balance
    der(HU_l*u_l+HU_v*u_v)=Vdot_in_l*rho_l_in/MM_l_in*h_l_in-Vdot_out_l*rho_l/MM_l*h_l-Vdot_out_v*rho_v/MM_v*h_v+Qdot_in-Qdot_out
    "Wärmeaufnahme des Stahls fehlt noch";
    T_l=T_v;
    T_system=T_l;

/* heat transfer */
Qdot_in=heatTransferModel.Q_kon;

    p_out_l_port=p_system+rho_l*g*height_l;

    //p_system=p_out_v_port; //Ist für Schaltung mit Kolonne überflüssig

    //p_system=p_in;
        V_l=HU_l*MM_l/rho_l;
    V_v=HU_v*MM_v/rho_v;

    Vdot_out_v=(p_system-p_out_v_port)/1e3;

assert(p_system>p_out_v_port, "p<p_out; pressure at outlet higher than inside the system");
assert((HU_l*MM_l/rho_l/Vges)<0.95, "eps_liq>0.95; system ist flooded");

  annotation (Diagram(graphics), Icon(graphics={
        Ellipse(
          extent={{-62,18},{38,-64}},
          lineColor={0,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={0,128,255}),
        Rectangle(
          extent={{-62,70},{38,16}},
          lineColor={0,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={170,213,255}),
        Ellipse(
          extent={{-62,64},{38,78}},
          lineColor={0,0,0},
          fillPattern=FillPattern.CrossDiag,
          fillColor={170,213,255}),
        Rectangle(
          extent={{-62,24},{38,-30}},
          lineColor={0,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={0,128,255}),
        Rectangle(
          extent={{-68,-24},{38,-30}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={255,0,0}),
        Rectangle(
          extent={{-62,24},{44,18}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={255,0,0}),
        Rectangle(
          extent={{-62,18},{-56,-24}},
          lineColor={0,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={255,0,0}),
        Rectangle(
          extent={{-61,21},{-57,-27}},
          fillPattern=FillPattern.Solid,
          fillColor={255,0,0},
          pattern=LinePattern.None),
        Rectangle(
          extent={{-46,18},{-40,-24}},
          lineColor={0,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={255,0,0}),
        Rectangle(
          extent={{-45,21},{-41,-27}},
          fillPattern=FillPattern.Solid,
          fillColor={255,0,0},
          pattern=LinePattern.None),
        Rectangle(
          extent={{-30,18},{-24,-24}},
          lineColor={0,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={255,0,0}),
        Rectangle(
          extent={{-29,21},{-25,-27}},
          fillPattern=FillPattern.Solid,
          fillColor={255,0,0},
          pattern=LinePattern.None),
        Rectangle(
          extent={{-14,18},{-8,-24}},
          lineColor={0,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={255,0,0}),
        Rectangle(
          extent={{-13,21},{-9,-27}},
          fillPattern=FillPattern.Solid,
          fillColor={255,0,0},
          pattern=LinePattern.None),
        Rectangle(
          extent={{2,18},{8,-24}},
          lineColor={0,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={255,0,0}),
        Rectangle(
          extent={{3,21},{7,-27}},
          fillPattern=FillPattern.Solid,
          fillColor={255,0,0},
          pattern=LinePattern.None),
        Rectangle(
          extent={{18,18},{24,-24}},
          lineColor={0,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={255,0,0}),
        Rectangle(
          extent={{19,21},{23,-27}},
          fillPattern=FillPattern.Solid,
          fillColor={255,0,0},
          pattern=LinePattern.None),
        Rectangle(
          extent={{32,18},{38,-24}},
          lineColor={0,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={255,0,0}),
        Rectangle(
          extent={{33,21},{37,-27}},
          fillPattern=FillPattern.Solid,
          fillColor={255,0,0},
          pattern=LinePattern.None)}),
    Documentation(info="<html>
<p><h4>Verdampfer</h4></p>
<p>Der Verdampfer wird so modelliert, dass er einen Eingang hat, in den die komplette Fl&uuml;ssigkeit eintritt und zwei Austritte. &Uuml;ber den ersten wird die verdampfte Fl&uuml;ssigkeit wieder in die Kolonne gef&uuml;hrt und &uuml;ber den zweiten wird das Sumpfprodukt, also die nicht verdampfte Fl&uuml;ssigkeit abgezogen.</p>
<p>Dieser Ausfluss wird geregelt, sodass im Verdampfer immer ausreichend Fl&uuml;ssigkeit vorhanden ist.</p>
<p>Der Inhalt des Verdampfers besteht aus einer Gasphase und einer Fl&uuml;ssigphase, die im thermodynamischen Gleichgewicht stehen.</p>
<p>Es wird bei der Modellierung angenommen, dass der Siedepunkt des Leichtsieders, in diesem Fall Ethanol, schon erreicht ist, da sonst keine Gasphase existieren w&uuml;rde, und somit die Gleichung f&uuml;r das Phasengleichgewicht nicht anwendbar w&auml;re. Erweiterungen diesbez&uuml;glich werden in Kapitel \\ref{subsec:Erw.2} vorgestellt. Es werden also vorerst keine Anfahrvorg&auml;nge aus dem kalten und leeren Zustand betrachtet. Der Sumpf ist von Beginn an bis zur Sollgrenze mit Fl&uuml;ssigkeit gef&uuml;llt, die bereits eine Temperatur von etwa 95 &deg; C hat.</p>
<p>Zur Berechnung der Austrittsgr&ouml;&szlig;en, die in diesem Modell die Unbekannten darstellen, stehen die Eintrittsgr&ouml;&szlig;en, sowie einige Gr&ouml;&szlig;en, die vom Benutzer gew&auml;hlt werden, zur Verf&uuml;gung. Es werden dabei folgende Gr&ouml;&szlig;en berechnet:\\par\\medskip\\noindent</p>
<p>Gasaustritt</p>
<p><ul>
<li> Volumenstrom <i>Vdot_V</i></li>
<li> Zusammensetzung <i>x_V</i></li>
<li>Temperatur <i>T_out_V</i></li>
<li> Druck p_out_V</li>
<li> Enthalpie <i>h_out_V</i></li>
</ul></p>
<p><br/>Fl&uuml;ssigkeitsaustritt</p>
<p><ul>
<li>Volumenstrom <i>Vdot_L</i></li>
<li>Zusammensetzung<i> x_L</i></li>
<li>Temperatur <i>T_out_L</i></li>
<li>Druck <i>p_out_L</i></li>
<li>Enthalpie <i>h_out_L</i></li>
</ul></p>
<p>Die gleichen Gr&ouml;&szlig;en, die unbekannt f&uuml;r die Austritte sind, sind vorgegeben aus dem untersten Boden der Kolonne f&uuml;r den Eintritt in den Verdampfer. Die zugef&uuml;hrte W&auml;rme <i>Qdot_in</i>, sowie die sich aus den Ma&szlig;en des Beh&auml;lters ergebenden Geometrievariablen, also die Fl&auml;che <i>A </i>und die H&ouml;he <i>h </i>werden vorgegeben. Au&szlig;erdem wird wird der gew&uuml;nschte F&uuml;llstand, &uuml;ber den der Ausfluss der Fl&uuml;ssigkeit geregelt wird <i>h_w</i>, festgelegt. </p>
<p>Eine Verlustw&auml;rme nach au&szlig;en zur Umgebung wird vernachl&auml;ssigt. Weitere Annahmen sind, dass der Eintrittsdruck der Fl&uuml;ssigkeit dem Systemdruck im Verdampfer entspricht und somit auch dem Austrittsdruck des Gases.</p>
</html>"));
end Sump_Steam_1_2;
