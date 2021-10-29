within ThermalSeparation.Components.Condenser;
model TotalCondenser_PT1
extends ThermalSeparation.Icons.Color.Condenser;
parameter Real timeConstant = 100;
    outer ThermalSeparation.SystemTS systemTS;
      parameter Modelica.SIunits.Temperature T_ref=systemTS.T_ref
    "reference temperature"                                                 annotation (Dialog(tab="Advanced"));
replaceable package MediumVapour =
      ThermalSeparation.Media.C2H5OH_Water_Vap                                                  constrainedby ThermalSeparation.Media.BaseMediumVapour
                                                                                                        annotation(choicesAllMatching);

    MediumVapour.BaseProperties mediumVapourIn(T0 = T_ref, c=x_v*rho_v/MM_v,p=vapourIn.p, T=T_v, x=x_v,  x_star=x_v);

  replaceable package MediumLiquid =
  ThermalSeparation.Media.C2H5OH_Water_Liq            constrainedby ThermalSeparation.Media.BaseMediumLiquid
                                                                                                        annotation(choicesAllMatching);

    MediumLiquid.BaseProperties mediumLiquid(T0 = T_ref, p=p_cond, T=T_l, x=x_l);
    MediumLiquid.BaseProperties mediumLiquid2(T0 = T_ref, p=p_cond, T=T_cond, x=x_l);

parameter Integer nS=2;
final parameter Integer nL=nSL-nS;
final parameter Integer nV=nSV-nS;
  final parameter Integer nSL = MediumLiquid.nSubstance;
  final parameter Integer nSV = MediumVapour.nSubstance;

parameter Boolean inertVapour[nSV] = fill(false,nSV)
    "true for each component which is inert";
parameter Boolean inertLiquid[nSL] = fill(false,nSL)
    "true for each component which is inert";

//inflow

    // SI.Pressure p_v;//=2e5;
     Modelica.SIunits.MoleFraction x_v[nSV];

     Modelica.SIunits.Temperature T_v;

     ThermalSeparation.Units.MolarEnthalpy h_v;
     ThermalSeparation.Units.MolarEnthalpy u_v;
     Modelica.SIunits.Density rho_v;
     Modelica.SIunits.MolarMass MM_v;

parameter Modelica.SIunits.Pressure p_cond=1e5 "Kondensationsdruck";
     Modelica.SIunits.MoleFraction x_l[nSL];
    // SI.MoleFraction x_v_out[nSV];
     Modelica.SIunits.Temperature T_l;

     ThermalSeparation.Units.MolarEnthalpy h_l;
     ThermalSeparation.Units.MolarEnthalpy u_l;
     Modelica.SIunits.Density rho_l;
     Modelica.SIunits.MolarMass MM_l;
     //SI.Pressure p_sat[nS];

//Druck in der obersten Stufe
Modelica.SIunits.Pressure p_1;

//SI.Temperature T_set; // feste Temperatur für das Kondensat
parameter Modelica.SIunits.Temperature T_set=303.15;
Modelica.SIunits.HeatFlowRate Q_kond;

      ThermalSeparation.Interfaces.LiquidPortOut liquidOut(redeclare package
      Medium = MediumLiquid)
        annotation (Placement(transformation(extent={{-6,-88},{14,-68}}),
        iconTransformation(extent={{-26,-108},{14,-68}})));
   ThermalSeparation.Interfaces.GasPortIn vapourIn(redeclare package Medium =
        MediumVapour)
     annotation (Placement(transformation(extent={{-6,88},{14,108}}),
        iconTransformation(extent={{-26,68},{14,108}})));
/*** für die Bestimmung der Kondensationstemperatur ***/
  parameter
    ThermalSeparation.Components.Condenser.Enumerations.OutletTempOption
    outletTempOption=
         Enumerations.OutletTempOption.T_set
    "different options for liquid outlet temperature"
    annotation(Dialog(Evaluate=true));

parameter Integer mapping[nS,2]= {{1,1},{2,2}}
    "parameter to map the different medium vectors one to another";
parameter Real factor_K[nSL] = fill(1,nSL);
 MediumLiquid.ActivityCoefficient activityCoeff(T=T_cond,x_l=x_l);
 MediumLiquid.FugacityCoefficient fugacityCoeffSat(T=T_cond, p=p_cond, p_sat=mediumLiquid2.p_sat);
 Modelica.SIunits.MoleFraction x_v_star[nSV];
 Modelica.SIunits.Temperature T_cond(start=63 + 273);
 parameter Modelica.SIunits.Temperature delta_T_sub=5 "subcooling";

  Modelica.Blocks.Continuous.FirstOrder firstOrder(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=65 + 273,
    T=timeConstant)
    annotation (Placement(transformation(extent={{-58,46},{-38,66}})));

  Modelica.Blocks.Continuous.FirstOrder firstOrder1(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=65 + 273,
    T=timeConstant)
    annotation (Placement(transformation(extent={{62,46},{82,66}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder2(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=65 + 273,
    T=timeConstant)
    annotation (Placement(transformation(extent={{-20,46},{0,66}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder3(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{-58,6},{-38,26}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder4(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{56,6},{76,26}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder5(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{-20,6},{0,26}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder6(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{-58,-26},{-38,-6}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder7(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{56,-26},{76,-6}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder8(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{-20,-26},{0,-6}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder9(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=1,
    T=timeConstant)
    annotation (Placement(transformation(extent={{-58,-56},{-38,-36}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder10(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=1,
    T=timeConstant)
    annotation (Placement(transformation(extent={{56,-56},{76,-36}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder11(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=1,
    T=timeConstant)
    annotation (Placement(transformation(extent={{-20,-56},{0,-36}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder12(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{-58,-84},{-38,-64}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder13(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{54,-84},{74,-64}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder14(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{-22,-84},{-2,-64}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder15(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{20,6},{40,26}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder16(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{20,-26},{40,-6}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder17(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=1,
    T=timeConstant)
    annotation (Placement(transformation(extent={{20,-56},{40,-36}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder18(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=0,
    T=timeConstant)
    annotation (Placement(transformation(extent={{18,-84},{38,-64}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder19(
    initType=Modelica.Blocks.Types.Init.InitialState,
    k=1,
    y_start=65 + 273,
    T=timeConstant)
    annotation (Placement(transformation(extent={{22,46},{42,66}})));
equation
  firstOrder.u = T_v;
firstOrder1.y = T_l;
firstOrder3.u = x_v[1];
firstOrder6.u = x_v[2];
firstOrder9.u = x_v[3];
firstOrder12.u = x_v[4];
firstOrder4.u = x_l[1];
firstOrder7.u = x_l[2];
firstOrder10.u = x_l[3];
firstOrder13.u = x_l[4];
   //Liquid inflow
      vapourIn.p = p_1;

       inStream(vapourIn.h_outflow) = h_v;
       inStream(vapourIn.x_outflow) = x_v;

   //Liquid outflow

       liquidOut.h_outflow = h_l;
       liquidOut.x_outflow = x_l;

    vapourIn.h_outflow = h_l;
       vapourIn.x_outflow = x_l;
//Medienübergabe

   MM_v = mediumVapourIn.MM;
   rho_v = mediumVapourIn.d;
   h_v = mediumVapourIn.h;
   u_v = mediumVapourIn.u;
  h_l =mediumLiquid.h;
  u_l =mediumLiquid.u;
  rho_l = mediumLiquid.d;
  MM_l = mediumLiquid.MM;

vapourIn.Ndot+liquidOut.Ndot=0;

p_1=p_cond;//zeta*rho_v/2*(w_v)^2+p_kond;//Druckverlust

//Energy balance
  Q_kond=vapourIn.Ndot*h_v +liquidOut.Ndot*h_l;
  //Q_kond=Vdot_out_CW*rho_CW_out/MM_CW_out*h_CW_out-Vdot_in_CW*rho_CW_in/MM_CW_in*h_CW_in;

//  Vdot_out_CW*rho_CW_out/MM_CW_out=Vdot_in_CW*rho_CW_in/MM_CW_in;

/*** Bestimmung der Kondensationstemp. wie bei Forner S. 48 ff ***/
for i in 1:nS loop
     x_v_star[mapping[i,1]] * p_cond *1 = factor_K[mapping[i,2]]* x_l[mapping[i,2]] *activityCoeff.gamma[mapping[i,2]] *fugacityCoeffSat.phi_sat[mapping[i,2]] * mediumLiquid2.p_sat[mapping[i,2]];
end for;
  sum(x_v_star)=1;
  connect(firstOrder.y, firstOrder2.u) annotation (Line(
      points={{-37,56},{-22,56}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder3.y, firstOrder5.u)
                                       annotation (Line(
      points={{-37,16},{-22,16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder6.y, firstOrder8.u)
                                       annotation (Line(
      points={{-37,-16},{-22,-16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder9.y, firstOrder11.u)
                                       annotation (Line(
      points={{-37,-46},{-22,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder12.y, firstOrder14.u)
                                       annotation (Line(
      points={{-37,-74},{-24,-74}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder5.y, firstOrder15.u)
                                        annotation (Line(
      points={{1,16},{18,16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder8.y, firstOrder16.u)
                                        annotation (Line(
      points={{1,-16},{18,-16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder11.y,firstOrder17. u)
                                        annotation (Line(
      points={{1,-46},{18,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder15.y, firstOrder4.u) annotation (Line(
      points={{41,16},{54,16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder16.y, firstOrder7.u) annotation (Line(
      points={{41,-16},{48,-16},{48,-16},{54,-16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder17.y, firstOrder10.u) annotation (Line(
      points={{41,-46},{46,-46},{46,-46},{54,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder18.y, firstOrder13.u) annotation (Line(
      points={{39,-74},{52,-74}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder2.y, firstOrder19.u) annotation (Line(
      points={{1,56},{20,56}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder19.y, firstOrder1.u) annotation (Line(
      points={{43,56},{60,56}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder14.y, firstOrder18.u) annotation (Line(
      points={{-1,-74},{16,-74}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics), Icon(coordinateSystem(preserveAspectRatio=
            false, extent={{-100,-100},{100,100}}),
                                      graphics={
        Text(
          extent={{-43,1.45455},{43,25.4545}},
          lineColor={0,0,0},
          origin={-6.5455,67},
          rotation=180,
          textString="%name")}),
    Documentation(info="<html>
<p><h4><font color=\"#008000\">Kondensator</font></h4></p>
<p>Der Kondensator hat einen Eingang in den das dampff&ouml;rmige Medium aus der obersten Stufe der Kolonne eintritt sowie einen Ausgang, durch den der komplett kondensierte Dampf als Fl&uuml;ssigkeit teilweise wieder in die Kolonne zur&uuml;ckgef&uuml;hrt wird.</p>
<p>Zwischen Austritt des Kondensators und der Kolonne wird ein Splitter gesetzt, &uuml;ber den das R&uuml;cklaufverh&auml;ltnis festgelegt wird.</p>
<p>Das Modell wird als einfacher Bilanzraum mit zwei Volumenstr&ouml;men und einem aus dem Kondensator austretenden W&auml;rmestrom betrachtet.</p>
<p>Es wird bei der Modellentwicklung Totalkondensation angenommen. </p>
<p>Das hei&szlig;t, dass der eintretende reine Dampfstrom den Kondensator als Fl&uuml;ssigkeit verl&auml;sst.</p>
<p>Au&szlig;er dem Druck sind alle Gr&ouml;&szlig;en des eintretenden Volumenstroms aus der obersten Stufe der Kolonne festgelegt. Es sind somit eintretender Volumenstrom <i>Vdot_in_V</i>, Zusammensetzung <i>y_in</i>, Temperatur <i>T_in_V</i> und spezifische Enthalpie <i>h_in_V</i> bekannt.</p>
<p>Der Druck wird f&uuml;r das Gesamtmodell im Kondensatormodell vorgegeben. Die weiteren Dr&uuml;cke der einzelnen B&ouml;den werden durch Druckverlustgleichungen berechnet.</p>
<p>Der Druckverlust vom Kondensator zur obersten Stufe der Kolonne wird hierbei vernachl&auml;ssigt.</p>
<p>F&uuml;r den austretenden fl&uuml;ssigen Strom m&uuml;ssen also die Gr&ouml;&szlig;en f&uuml;r folgende Variablen bestimmt werden:</p>
<p><ul>
<li> Volumenstrom <i>Vdot_out_L</i></li>
<li> Zusammensetzung <i>x_out_L</i></li>
<li> Temperatur <i>T_out_L</i></li>
<li> Druck <i>p_out_L</i></li>
<li> Enthalpie <i>h_out_L</i></li>
</ul></p>
</html>"));
end TotalCondenser_PT1;
