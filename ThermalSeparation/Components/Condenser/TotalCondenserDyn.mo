within ThermalSeparation.Components.Condenser;
model TotalCondenserDyn
extends ThermalSeparation.Icons.Color.Condenser;

    outer ThermalSeparation.SystemTS systemTS;
      parameter Modelica.SIunits.Temperature T_ref=systemTS.T_ref
    "reference temperature"                                                 annotation (Dialog(tab="Advanced"));
replaceable package MediumVapour =
      ThermalSeparation.Media.C2H5OH_Water_Vap                                                  constrainedby
    ThermalSeparation.Media.BaseMediumVapour                                                            annotation(choicesAllMatching);

    MediumVapour.BaseProperties mediumVapourIn(T0 = T_ref, c=c_v,p=vapourIn.p, T=T_v, x=x_v,  x_star=x_v);

  replaceable package MediumLiquid =
  ThermalSeparation.Media.C2H5OH_Water_Liq            constrainedby
    ThermalSeparation.Media.BaseMediumLiquid                                                            annotation(choicesAllMatching);

    MediumLiquid.BaseProperties mediumLiquid(T0 = T_ref, p=p_cond, T=T_l, x=x_l,h=h_l);
//    MediumLiquid.BaseProperties mediumLiquid2(T0 = T_ref, p=p_cond, T=T_cond, x=x_l);

parameter Integer nS=2;
final parameter Integer nL=nSL-nS;
final parameter Integer nV=nSV-nS;
  final parameter Integer nSL = MediumLiquid.nSubstance;
  final parameter Integer nSV = MediumVapour.nSubstance;

parameter Boolean inertVapour[nSV] = fill(false,nSV)
    "true for each component which is inert";
parameter Boolean inertLiquid[nSL] = fill(false,nSL)
    "true for each component which is inert";

     Modelica.SIunits.MoleFraction x_v[nSV];

     Modelica.SIunits.Temperature T_v;

     Modelica.SIunits.Concentration c_v[nSV]=x_v*rho_v/MM_v;
     ThermalSeparation.Units.MolarEnthalpy h_v= mediumVapourIn.h;
     ThermalSeparation.Units.MolarEnthalpy u_v= mediumVapourIn.u;
     Modelica.SIunits.Density rho_v= mediumVapourIn.d;
     Modelica.SIunits.MolarMass MM_v= mediumVapourIn.MM;

parameter Modelica.SIunits.Pressure p_cond=1e5 "Kondensationsdruck";
     Modelica.SIunits.MoleFraction x_l[nSL];
    // SI.MoleFraction x_v_out[nSV];
     Modelica.SIunits.Temperature T_l;
     ThermalSeparation.Units.MolarEnthalpy h_l;//=mediumLiquid.h;

//Druck in der obersten Stufe
Modelica.SIunits.Pressure p_1;

parameter Modelica.SIunits.Temperature T_set=303.15;
Modelica.SIunits.HeatFlowRate Q_kond;

      ThermalSeparation.Interfaces.LiquidPortOut liquidOut(redeclare package
      Medium = MediumLiquid)
        annotation (Placement(transformation(extent={{-6,-86},{14,-66}}),
        iconTransformation(extent={{-26,-108},{14,-66}})));
   ThermalSeparation.Interfaces.GasPortIn vapourIn(redeclare package Medium =
        MediumVapour)
     annotation (Placement(transformation(extent={{-6,88},{14,108}}),
        iconTransformation(extent={{-26,68},{14,108}})));
/*** fr die Bestimmung der Kondensationstemperatur ***/
  parameter
    ThermalSeparation.Components.Condenser.Enumerations.OutletTempOption
    outletTempOption=
         Enumerations.OutletTempOption.T_set
    "different options for liquid outlet temperature"
    annotation(Dialog(Evaluate=true));

parameter Integer mapping[nS,2]= {{1,1},{2,2}}
    "parameter to map the different medium vectors one to another";
parameter Real factor_K[nSL] = fill(1,nSL);
// MediumLiquid.ActivityCoefficient activityCoeff(T=T_cond,x_l=x_l);
// MediumLiquid.FugacityCoefficient fugacityCoeffSat(T=T_cond, p=p_cond, p_sat=mediumLiquid2.p_sat);
// Modelica.SIunits.MoleFraction x_v_star[nSV];
// Modelica.SIunits.Temperature T_cond(start=63 + 273);
 parameter Modelica.SIunits.Temperature delta_T_sub=5 "subcooling";

 parameter SI.Volume V_liq = 0.1 "constant liquid holdup in condenser";
 parameter SI.MoleFraction x_l_start[nSL];

equation
   //inflow
      vapourIn.p = p_1;
       inStream(vapourIn.h_outflow) = h_v;
       inStream(vapourIn.x_outflow) = x_v;
vapourIn.x_outflow =x_l;
vapourIn.h_outflow =h_l;
   //outflow

      // liquidOut.p = p_1;
       liquidOut.h_outflow = h_l;
       liquidOut.x_outflow = x_l;

  for i in 1:nS loop
     V_liq*der(mediumLiquid.c[i]) = vapourIn.Ndot*x_v[i] +liquidOut.Ndot * x_l[i];
//   inStream(vapourIn.x_outflow[i]) = (liquidOut.x_outflow[i]);
  end for;
 vapourIn.Ndot + liquidOut.Ndot=0;

  //sum(x_v_out)=1;
    if outletTempOption == Enumerations.OutletTempOption.T_in then
    T_l=mediumVapourIn.T;
    elseif outletTempOption ==  Enumerations.OutletTempOption.T_set then
      T_l = T_set;
   elseif outletTempOption ==  Enumerations.OutletTempOption.T_subcool then
     T_l=300;//T_cond - delta_T_sub;
    end if;

//Druckberechnung fr die oberste Stufe

p_1=p_cond;//zeta*rho_v/2*(w_v)^2+p_kond;//Druckverlust

//Energy balance
  Q_kond= vapourIn.Ndot*inStream(vapourIn.h_outflow) + liquidOut.Ndot*(liquidOut.h_outflow);

/*** Bestimmung der Kondensationstemp. wie bei Forner S. 48 ff ***/
for i in 1:nS loop
 //    x_v_star[mapping[i,1]] * p_cond *1 = factor_K[mapping[i,2]]* x_l[mapping[i,2]] *activityCoeff.gamma[mapping[i,2]] *fugacityCoeffSat.phi_sat[mapping[i,2]] * mediumLiquid2.p_sat[mapping[i,2]];
end for;
 // sum(x_v_star)=1;

initial equation
  x_l=x_l_start;
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
end TotalCondenserDyn;
