within ThermalSeparation.Components.ColumnSimp;
model OneStageSimp

extends ThermalSeparation.Components.Condenser.BaseCondenser(connectCoolant=true);

extends ThermalSeparation.Icons.Color.Condenser;

MediumVapour.BaseProperties mediumVapour(T=T_v, x=x_v, p=p, x_star=x_v, c=c_v);
MediumLiquid.BaseProperties mediumLiquidIn(T=T_l_in,x=x_l_in, p=p,h=h_l_in);

parameter Boolean inert_Vapour[nSV] = fill(false,nSV)
    "true for inert components in vapour phase" annotation(Dialog(group="Shell side parameters"));
parameter Boolean inert_Liquid[nSL] = fill(false,nSL)
    "true for inert components in vapour phase" annotation(Dialog(group="Shell side parameters"));

  Real K[nS] = thermoEquilibrium.K "equilibrium constant";

  parameter Modelica.Units.SI.HeatFlowRate Q_kon_fixed=1e5;
  parameter Real beta_l=7e2 "liquid mass transfer coefficient" annotation(Dialog(group="Shell side parameters"));
  parameter Real beta_v=7e2 "vapour mass transfer coefficient" annotation(Dialog(group="Shell side parameters"));
  parameter Real alpha_l=5e5 "liquid heat transfer coefficient" annotation(Dialog(group="Shell side parameters"));
  parameter Real alpha_v=5e5 "vapour heat transfer coefficient" annotation(Dialog(group="Shell side parameters"));

  Modelica.Units.SI.MolarFlowRate ndot_l_interface[nS];
  Modelica.Units.SI.MolarFlowRate ndot_v_interface[nS];
  Modelica.Units.SI.MolarFlowRate ndot_from_l[nSL];
  Modelica.Units.SI.MolarFlowRate ndot_from_v[nSV];

  Modelica.Units.SI.Temperature T_l_in;
  Modelica.Units.SI.Temperature T_v(stateSelect=StateSelect.default);
  Modelica.Units.SI.Temperature Tstar;

  Modelica.Units.SI.MoleFraction x_l_in[nSL];
  Modelica.Units.SI.MoleFraction x_v[nSV];
  Modelica.Units.SI.MoleFraction x_l_star[nSL];
  Modelica.Units.SI.MoleFraction x_v_star[nSV];
  Modelica.Units.SI.MoleFraction x_from_l[nSL];
  Modelica.Units.SI.MoleFraction x_to_l[nSL];
  Modelica.Units.SI.MoleFraction x_from_v[nSV];
  Modelica.Units.SI.MoleFraction x_to_v[nSV];

  Modelica.Units.SI.MolarMass MM_l_in=mediumLiquidIn.MM;
  Modelica.Units.SI.MolarMass MM_v=mediumVapour.MM;

  Modelica.Units.SI.Density rho_l_in=mediumLiquidIn.d;
  Modelica.Units.SI.Density rho_v=mediumVapour.d;

  Modelica.Units.SI.Concentration c_l_in[nSL];
  Modelica.Units.SI.Concentration c_v[nSV];

  ThermalSeparation.Units.MolarEnthalpy h_l_in;//mediumLiquidIn.h;
  ThermalSeparation.Units.MolarEnthalpy h_v = mediumVapour.h;
  ThermalSeparation.Units.MolarEnthalpy u_v = mediumVapour.u;
  ThermalSeparation.Units.MolarEnthalpy h_from_v=mediumVapour.h;//enthalpy_VapToPB.h;
  ThermalSeparation.Units.MolarEnthalpy h_from_l=mediumLiquid.h;//enthalpy_LiqToPB.h;
  ThermalSeparation.Units.MolarEnthalpy h_to_v=mediumVapour.h;//enthalpy_VapToPB.h;
  ThermalSeparation.Units.MolarEnthalpy h_to_l=mediumLiquid.h;//enthalpy_LiqToPB.h;

  Modelica.Units.SI.HeatFlowRate qdot_v_interface;
  Modelica.Units.SI.HeatFlowRate qdot_l_interface;
  Modelica.Units.SI.HeatFlowRate qdot_v_conv;
  Modelica.Units.SI.HeatFlowRate qdot_l_conv;
  Modelica.Units.SI.HeatFlowRate qdot_v_cond;
  Modelica.Units.SI.HeatFlowRate qdot_l_cond;

  Modelica.Units.SI.VolumeFlowRate vdot_v_out;
  Modelica.Units.SI.VolumeFlowRate vdot_l_in;

  SI.MoleFraction x_vap_liq[nS] "total molar fractions";
  Real n_tot[nS] "total molar holdup";

  Real eps_liq;
  Real eps_liq_init=0.25 "initial liquid volume fraction" annotation(Dialog(tab="Initialisation"));
  Real eps_vap;

  Real z;

/* medium properties */

/* geometry */

  parameter Modelica.Units.SI.Area A=0.5 "crossectional area" annotation(Dialog(tab="Geometry"));
  parameter Modelica.Units.SI.Height H=1.3 "height" annotation(Dialog(tab="Geometry"));

  Modelica.Units.SI.Height h_liq "liquid level inside condenser";
  parameter Modelica.Units.SI.Length h_w=H/10 "weir height" annotation(Dialog(tab="Geometry"));
  parameter Modelica.Units.SI.Length h_lw=A*0.7 "weir length" annotation(Dialog(tab="Geometry"));

MediumLiquid.FugacityCoefficient satFugacityLiq(T=Tstar, p=p, p_sat=p_sat);

   replaceable model ThermoEquilibrium =
     ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium constrainedby ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium
    "model for phase equilibrium" annotation (choicesAllMatching=true,
      Dialog(group="Shell side parameters"));

      ThermoEquilibrium thermoEquilibrium(nS=nS,mapping = mapping,
      redeclare replaceable package MediumVapour = MediumVapour, redeclare replaceable package MediumLiquid =
      MediumLiquid, p=p, T=Tstar, x_v=x_v_star, x_l=x_l_star, p_sat=p_sat,  v_v=MM_v./rho_v, x_vap_liq=x_vap_liq);

/* connectors */

     ThermalSeparation.Interfaces.GasPortIn gasPortIn(
       redeclare package Medium = MediumVapour)                                                                                                     annotation (Placement(transformation(extent={{-102,12},
            {-82,32}}),        iconTransformation(extent={{-102,12},{
            -82,32}})));
    ThermalSeparation.Interfaces.LiquidPortOut liquidPortOut(
      redeclare package Medium = MediumLiquid)
      annotation (Placement(transformation(extent={{70,-50},{90,-30}}),
          iconTransformation(extent={{70,-50},{90,-30}})));
    ThermalSeparation.Interfaces.GasPortOut gasPortOut(
      redeclare package Medium = MediumVapour)
      annotation (Placement(transformation(extent={{70,30},{90,50}}),
          iconTransformation(extent={{70,30},{90,50}})));
  ThermalSeparation.Interfaces.HeatPort heatPort(Qdot=-Q_kon)
    annotation (Placement(transformation(extent={{52,0},{72,20}}, rotation=0),
        iconTransformation(extent={{-10,70},{10,90}})));

    ThermalSeparation.Interfaces.LiquidPortIn liquidPortIn(redeclare package Medium =
               MediumLiquid)
      annotation (Placement(transformation(extent={{-100,-100},{-80,-80}}),
          iconTransformation(extent={{-80,-84},{-60,-64}})));
equation
/** Connector variables **/

inStream(gasPortIn.x_outflow)=x_v_in;
inStream(gasPortIn.h_outflow)=h_v_in;
gasPortIn.Ndot=vdot_v_in*rho_v_in/MM_v_in;
gasPortIn.p=p_in;
gasPortIn.x_outflow=x_v;
gasPortIn.h_outflow=h_v;

liquidPortOut.Ndot=-vdot_l_out*rho_l/MM_l;
liquidPortOut.p=p;
liquidPortOut.x_outflow=x_l;
liquidPortOut.h_outflow=h_l;

gasPortOut.Ndot=-vdot_v_out*rho_v/MM_v;
gasPortOut.p=p;
gasPortOut.x_outflow=x_v;
gasPortOut.h_outflow=h_v;
mdot_coolant_in=1; //dummy
Q_kon=Q_kon_fixed;

liquidPortIn.h_outflow = h_l;
liquidPortIn.x_outflow = x_l;
inStream(liquidPortIn.h_outflow) = h_l_in;
inStream(liquidPortIn.x_outflow) = x_l_in;
liquidPortIn.p = p;
liquidPortIn.Ndot = vdot_l_in*rho_l_in/MM_l_in;

  /*pressure*/
  vdot_v_out/vdot_nom=(p_in-p)/p_nom;

/* geometry */
   eps_liq=A*h_liq/(A*H);

   //der(eps_liq)+der(eps_vap)=0;
   eps_liq + eps_vap=1;

   vdot_l_out=1.848*h_lw*(abs(h_liq-h_w))^1.5; //francis weir formula

/* correlation between mole fraction and concentration */
  for i in 1:nSL loop
    c_l[i]=x_l[i]*rho_l/MM_l;
    c_l_in[i]=x_l_in[i]*rho_l_in/MM_l_in;
  end for;

  for i in 1:nSV loop
    c_v[i]=x_v[i]*rho_v/MM_v;
    c_v_in[i]=x_v_in[i]*rho_v_in/MM_v_in;
  end for;

/* mole balance */

  for i in 1:nSV loop
  /* mole balance vapour phase (by components) */
  A*H*der(eps_vap*c_v[i])=vdot_v_in*c_v_in[i]-vdot_v_out*c_v[i]+ndot_v_interface[i];

  /* mole balance for inert substances */
  if inert_Vapour[i] then
     A*H*der(eps_vap*c_v[i])=vdot_v_in*c_v_in[i]-vdot_v_out*c_v[i];
  end if;

  end for;
  for i in 1:nSL loop
  /* mole balance liquid phase (by components)*/
  A*H*der(eps_liq*c_l[i])=vdot_l_in*c_l_in[i]-vdot_l_out*c_l[i]+ndot_l_interface[i];

    /* mole balance for inert substances */
  if inert_Liquid[i] then
     A*H*der(eps_liq*c_l[i])=vdot_l_in*c_l_in[i]-vdot_l_out*c_l[i];
  end if;
  end for;

  /* total mole balance liquid*/
  A*H*der(eps_liq*rho_l/MM_l)=vdot_l_in*rho_l_in/MM_l_in-vdot_l_out*rho_l/MM_l+sum(ndot_l_interface[:]);
  /* total mole balance vapour*/
  A*H*der(eps_vap*rho_v/MM_v)=vdot_v_in*rho_v_in/MM_v_in-vdot_v_out*rho_v/MM_v+sum(ndot_v_interface[:]);

  /* mass transfer at interface*/
  for i in 1:nS-1 loop
    ndot_l_interface[mapping[i,2]]=-beta_l*(x_l[mapping[i,2]]-x_l_star[mapping[i,2]]);
    ndot_v_interface[mapping[i,1]]=-beta_v*(x_v[mapping[i,1]]-x_v_star[mapping[i,1]]);
  end for;

  for j in 1:nS loop
        ndot_l_interface[mapping[j,2]]=-ndot_v_interface[mapping[j,1]];
        //ndot_v_interface[j]=-1e1*(x_v[j]-x_v_star[j]);
  end for;

/* energy balance */
  /* energy balance vapour phase */
  A*H*der(eps_vap*u_v*sum(c_v[:]))=vdot_v_in*h_v_in*sum(c_v_in[:])-vdot_v_out*h_v*sum(c_v[:])+qdot_v_interface; //h_evap included in h_v
  /* energy balance liquid phase */
  A*H*der(eps_liq*u_l*sum(c_l[:]))=vdot_l_in*h_l_in*sum(c_l_in[:])-vdot_l_out*h_l*sum(c_l[:])+qdot_l_interface-Q_kon;
//A*H*der(eps_vap*u_v*sum(c_v[:])+eps_liq*u_l*sum(c_l[:]))=vdot_v_in*h_v_in*sum(c_v_in[:])-vdot_v_out*h_v*sum(c_v[:])-vdot_l_out*h_l*sum(c_l[:])-Q_kon;

  /* energy transfer at interface */
  qdot_v_interface=qdot_v_conv+qdot_v_cond;
  qdot_l_interface=qdot_l_conv+qdot_l_cond;

  qdot_v_cond= alpha_v*(Tstar-T_v);
  qdot_l_cond= alpha_l*(Tstar-T_l);

  qdot_l_conv=-sum(ndot_from_l[:])*h_from_l+sum(ndot_from_v[:])*h_to_l;
  qdot_v_conv=-sum(ndot_from_v[:])*h_from_v+sum(ndot_from_l[:])*h_to_v;

  for i in 1:nS loop
  ndot_from_l[mapping[i,2]]=-1*min(0,ndot_l_interface[mapping[i,2]]);
  ndot_from_v[mapping[i,1]]=-1*min(0,ndot_v_interface[mapping[i,1]]);

  x_from_l[mapping[i,2]]=ndot_from_l[mapping[i,2]]/max(1e-5,sum(ndot_from_l[mapping[i,2]]));
  x_to_l[mapping[i,2]]=x_from_v[mapping[i,2]];
  x_from_v[mapping[i,1]]=ndot_from_v[mapping[i,1]]/max(1e-5,sum(ndot_from_v[:]));
  x_to_v[mapping[i,1]]=x_from_l[mapping[i,2]];
  end for;

/*energy balance at phase boundary*/
  -qdot_v_interface-qdot_l_interface=0;//+sum(ndot_l_interface[:]*h_evap[:])=0;

/* summation equation at the phase boundary */
  sum(x_l_star[:])=1;
  sum(x_v_star[:])=1;

/* thermodynamic equilibrium */
   for i in 1:nS loop
     x_v_star[mapping[i,1]]= K[i] *x_l_star[mapping[i,2]];
   end for;

/* control of output */
  z=vdot_v_out;

/* for calculation of external equilibrium */

    x_vap_liq[:] = n_tot[:]./sum(n_tot[:]);
    for i in 1:nS loop
       n_tot[i]= (c_v[mapping[i,1]].* (1-eps_liq) + c_l[mapping[i,2]].* eps_liq)*H*A;
    end for;

// assert(p>p_out, "p<p_out; pressure at outlet higher than inside the system");
// assert(vdot_v_out>0, "negative vapour outflow");
// assert(eps_liq<0.95, "eps_liq>0.95; system ist flooded");

initial equation

      sum(x_l[:])=1;
      sum(x_v[:])=1;

     //T_v =273+90;
     T_l =273+90;
     //T_v=T_l;
     //vdot_v_out=vdot_v_in;
     //p=1e5;
     for i in 2:nS loop

     ndot_v_interface[mapping[i,1]]=0;

     end for;

     eps_liq=eps_liq_init;
     //eps_vap=1-eps_liq_init;

  annotation (Diagram(graphics), Icon(graphics),
    Documentation(info="<html>
<p>Check sagt: 4 Unbekannte zuviel. Das Modell luft aber trotzdem. Grund: Zhlweise der Unbekannten in den Konnektoren. Beispielsweise sieht der Input Konnektor nicht alle Gren als bekannt an. Er will entweder eine Flow- oder eine Potentialvariable haben und liefert entsprechend die andere. Die korrekte Formulierung ist aber unintuitiv und mllt den Quellcode voll.</p> 
</html>"));
end OneStageSimp;
