within ThermalSeparation.Components.Reboiler;
model KettleReboilerEq "equilibrium model"

  extends ThermalSeparation.Icons.Color.Reboiler;

  outer ThermalSeparation.SystemTS systemTS;

  parameter ThermalSeparation.Components.Reboiler.InitOptionEq
    init_option =                     Enumerations.InitializationOption.init_x "initialization options"                        annotation(Dialog(tab="Initialization"),Evaluate=true);
   parameter Real eps_liq_init "initial value for eps_liq" annotation(Dialog(tab = "Initialization"));
  parameter Modelica.SIunits.Temperature T_init "initial value fr Temperature"
                                                                 annotation(Dialog(tab = "Initialization"));
  parameter Modelica.SIunits.Pressure p_init "initial value for pressure"
                                annotation(Dialog(tab = "Initialization"));
  parameter Real fixed_mol_init "molarity of solvent" annotation(Dialog(tab="Initialization"));
  parameter Modelica.SIunits.MoleFraction x_total_start[nS] "mole fraction of component in vapour and liquid" annotation(Dialog(tab = "Initialization"));

  Modelica.SIunits.HeatFlowRate Q_in "heat flow to evaporate liquid";

  Modelica.SIunits.Temperature T_l_in;

  Modelica.SIunits.MoleFraction x_l_in[nSL];

  Modelica.SIunits.Concentration c_l_in[nSL];
  Modelica.SIunits.VolumeFlowRate vdot_l_in;

  parameter Boolean init_standalone = true "changes number of initial values if used in standalone operation";

  Modelica.SIunits.Pressure p_sys(start=2e5);
  Modelica.SIunits.MolarFlowRate F_in_l;
  parameter Modelica.SIunits.Temperature T_ref=systemTS.T_ref "reference temperature"
                                                                        annotation (Dialog(tab="Advanced"));
  parameter Modelica.SIunits.HeatFlowRate Q_loss=0 "heat loss to ambience";
                                      //adiabatic system

  parameter Integer mapping[nS,2] = {{1,1},{2,2}} "parameter to map the different medium vectors one to another";

  parameter Integer nS=2 "number of omspecies which are equal in vapour and liquid phase";
  final parameter Integer nL=nSL-nS "number of additional substances which are only in liquid phase";
  final parameter Integer nV= nSV-nS "number of additional substances which are only in the vapour phase";
  final parameter Integer nSL = MediumLiquid.nSubstance;
  final parameter Integer nSV = MediumVapour.nSubstance;

  parameter Boolean inert_Liquid[nSL] = fill(false,nSL) "true for inert components in liquid phase";

  Real K[nS](start=fill(1,nS));

  Modelica.SIunits.MolarFlowRate F_out_l(start=1e-4);
  Modelica.SIunits.MolarFlowRate F_out_v(start=1e-4);

  Modelica.SIunits.AmountOfSubstance HU_l(start=1000, stateSelect=StateSelect.prefer);
  Modelica.SIunits.AmountOfSubstance HU_v(start=100);

  replaceable package MediumVapour = ThermalSeparation.Media.C2H5OH_Water_Vap
    constrainedby ThermalSeparation.Media.BaseMediumVapour "medium to be used in vapour phase"
                                                                         annotation(choicesAllMatching);

    MediumVapour.BaseProperties mediumVapour(
    T0=T_ref,
    p=p_sys,
    T=T,
    x=x_v,
    c=c_v, x_star=x_v);

  replaceable package MediumLiquid = ThermalSeparation.Media.C2H5OH_Water_Liq
    constrainedby ThermalSeparation.Media.BaseMediumLiquid "medium to be used in liquid phase"                                                         annotation(choicesAllMatching);

  MediumLiquid.BaseProperties mediumLiquid(
    T0=T_ref,
    p=p_sys,
    T=T,
    x=x_l,
    h=h_l);

  MediumLiquid.BaseProperties mediumLiquidIn(
    T0=T_ref,
    p=p_sys,
    T=T_l_in,
    x=x_l_in,
    h=h_l_in);

  replaceable model ThermoEquilibrium =
      ThermalSeparation.PhaseEquilibrium.H2O_CO2_MEA
   constrainedby ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium "model for phase equilibrium"  annotation (choicesAllMatching=true);

      ThermoEquilibrium thermoEquilibrium(nS=nS,mapping = mapping,
      redeclare replaceable package MediumVapour = MediumVapour, redeclare
      replaceable package                                                                      MediumLiquid =
      MediumLiquid, p=p_sys, T=T, x_v=x_v, x_l=x_l, p_sat=p_sat,  v_v=MM_v./rho_v,x_vap_liq=fill(1/nS,nS));

/* Reaction */
Modelica.SIunits.MolarFlowRate Ndot_l_transfer[nSL];

  replaceable model Reaction =
  ThermalSeparation.Reaction.NoReaction constrainedby
    ThermalSeparation.Reaction.BaseReaction                                                    "model for chemical reaction"                               annotation(choicesAllMatching=true);

  Reaction reaction(
   final n=1, c=c_l, V=V_liq,Ndot_l_transfer=Ndot_l_transfer,  gamma=fill(1,nSL),
     redeclare package MediumLiquid=MediumLiquid, propsLiq=mediumLiquid.properties);

  constant Real R=Modelica.Constants.R;

  package Col_Geometry =
  ThermalSeparation.Geometry;
  Col_Geometry.BasicGeometry geometry;

/* inner heat transfer resistance */

  replaceable model InnerHT =
      ThermalSeparation.HeatAndMassTransfer.HTResistance.NoHTResistance
    constrainedby
    ThermalSeparation.HeatAndMassTransfer.HTResistance.BaseHTResistance               "heat transfer mechanism between bulk and wall"                     annotation(choicesAllMatching=true);

  InnerHT innerHT(n=1,T={T},A={A_HT},Qdot={Q_in},p={p_sys});

  parameter Modelica.SIunits.Area A_HT "heat exchange area";

/* medium properties */
  Modelica.SIunits.Temperature T;

  Modelica.SIunits.MoleFraction x_l[nSL](start=fill(1/nSL,nSL));
  Modelica.SIunits.MoleFraction x_v[nSV](start=fill(1/nSV,nSV));

  ThermalSeparation.Units.MolarEnthalpy h_l_in;// = mediumLiquidIn.h;
  ThermalSeparation.Units.MolarEnthalpy h_v = mediumVapour.h;
  ThermalSeparation.Units.MolarEnthalpy h_l;// = mediumLiquid.h;
  ThermalSeparation.Units.MolarEnthalpy h_feed = 0; //simplification
  ThermalSeparation.Units.MolarEnthalpy u_l(stateSelect=StateSelect.prefer) = mediumLiquid.u;
  ThermalSeparation.Units.MolarEnthalpy u_v = mediumVapour.u;

  Modelica.SIunits.MolarMass MM_l=mediumLiquid.MM;
  Modelica.SIunits.MolarMass MM_v(start=0.031)=mediumVapour.MM;
  Modelica.SIunits.MolarMass MM_l_in=mediumLiquidIn.MM;

  Modelica.SIunits.Density rho_l=mediumLiquid.d;
  Modelica.SIunits.Density rho_v=mediumVapour.d;
  Modelica.SIunits.Density rho_l_in=mediumLiquidIn.d;

  Modelica.SIunits.Pressure p_sat[nSL];

  Modelica.SIunits.Volume V_vap;
  Modelica.SIunits.Volume V_liq;
  Modelica.SIunits.Volume V_abs;

  Modelica.SIunits.Concentration c_l[nSL];
  Real dummy_c_l(stateSelect=StateSelect.prefer)=c_l[2];
  Modelica.SIunits.Concentration c_v[nSV];

  Modelica.SIunits.VolumeFlowRate vdot_l;
  Modelica.SIunits.VolumeFlowRate vdot_v(start=50);

/* geometry */

  parameter Modelica.SIunits.Area A=0.5 "area of reboiler";
  parameter Modelica.SIunits.Height H=1.3 "height of reboiler";
  Modelica.SIunits.SpecificHeatCapacity cp_col=500;
                                                //stainless steel
  Modelica.SIunits.Mass m_col=100;
                     //mass of column segment

  Modelica.SIunits.Height h_liq;
 parameter  Modelica.SIunits.Length h_w=0.13;
 parameter Modelica.SIunits.Length h_lw=0.35;

  Real eps_liq "(liquid volume) / (absolute volume)";
  Real eps_vap(min = 0, max = 1) "(vapour volume) / (absolute volume)";
  //Real eps_inert "(inert volume) / (absolute volume)";
  Modelica.SIunits.Pressure p_bub "bubble point pressure";

  Real checkMoleBal;

  //for monitoring purpose

  Modelica.SIunits.MassFlowRate mdot_l = vdot_l * rho_l;
  Modelica.SIunits.MassFlowRate mdot_v = vdot_v * rho_v;
  Modelica.SIunits.MassFlowRate mdot_l_in = vdot_l_in * rho_l_in;
  Modelica.SIunits.MassFlowRate mdot_check = mdot_l_in - mdot_l - mdot_v;

  ThermalSeparation.Interfaces.GasPortOut gasPortOut(redeclare package Medium =
        MediumVapour) annotation (Placement(transformation(extent={{-40,82},{-20,
            102}}), iconTransformation(extent={{-60,62},{-20,102}})));
  ThermalSeparation.Interfaces.LiquidPortIn liquidPortIn(redeclare package Medium =
               MediumLiquid) annotation (Placement(transformation(extent={{28,
            82},{48,102}}), iconTransformation(extent={{8,62},{48,102}})));
  ThermalSeparation.Interfaces.LiquidPortOut liquidPortOut(redeclare package Medium =
               MediumLiquid) annotation (Placement(transformation(extent={{-100,
            -100},{-80,-80}}), iconTransformation(extent={{-28,-108},{12,-68}})));

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation (
      Placement(transformation(extent={{72,-10},{92,10}}),
        iconTransformation(extent={{72,-10},{92,10}})));
parameter Boolean init_const_vdot=false "initialise from constant vapour volume flow rate" annotation(Dialog(tab="Initialization",group="Init. from constant value"));
parameter Modelica.SIunits.VolumeFlowRate vdot_v_const=0.08 "volume flow rate" annotation(Dialog(tab="Initialization",group="Init. from constant value"));
   parameter Real omega_k=0.05 "large value if change between constant and variable shall be steep"
                                                                         annotation(Dialog(tab="Initialization",group="Init. from constant value"));
   parameter Real omega_time = 50 "inflexion point of tanh function" annotation(Dialog(tab="Initialization",group="Init. from constant value"));
   Real omega = 0.5*tanh(omega_k*(time-omega_time))+0.5;

equation
    /*Ports*/

    gasPortOut.p = p_sys;
    gasPortOut.h_outflow = h_v;
    gasPortOut.x_outflow = x_v;
    if Q_in==0 and vdot_l_in<=0.0003 then F_out_v=0; else gasPortOut.Ndot = -F_out_v; end if;
    //gasPortOut.Ndot = -F_out_v;

    liquidPortIn.p = p_sys;
    liquidPortIn.Ndot = F_in_l;
    inStream(liquidPortIn.h_outflow) = h_l_in;
    inStream(liquidPortIn.x_outflow) = x_l_in;
    liquidPortIn.x_outflow = x_l;
    liquidPortIn.h_outflow = h_l;

    liquidPortOut.Ndot = -F_out_l;
    liquidPortOut.h_outflow = h_l;
    liquidPortOut.x_outflow = x_l;
  //   liquidPortOut.p = p_sys;

    heatPort.Q_flow = Q_in;
    heatPort.T = innerHT.Twall[1];
    // vdot_v = 100;

  /* correlation between mole fraction and concentration */
  for i in 1:nSL loop
    c_l[i]=x_l[i]*rho_l/MM_l;
    c_l_in[i]=x_l_in[i]*rho_l_in/MM_l_in;
  end for;

  for i in 1:nSV loop
    c_v[i]=x_v[i]*rho_v/MM_v;
  end for;

  /* geometry */
  V_vap = eps_vap*V_abs;
  V_abs = A * H;
  V_liq = eps_liq*V_abs;

  eps_vap = (1 - eps_liq);

  h_liq = V_liq / A;

  p_bub = sum(x_l .* p_sat);

  F_out_l = noEvent(if h_liq<=h_w then 0 else rho_l/MM_l*1.848*h_lw*(abs(h_liq-h_w))^1.5); //francis weir formula

  /* mole balance */

  for i in 1:nS loop
   //der(V_liq*c_l[i]+V_vap*c_v[i]) = vdot_l_in * c_l_in[i] - vdot_l * c_l[i] - vdot_v * c_v[i] + reaction.Ndot[i];
      der(V_liq*c_l[mapping[i,2]]+V_vap*c_v[mapping[i,1]]) = liquidPortIn.Ndot * x_l_in[mapping[i,2]] + liquidPortOut.Ndot * x_l[mapping[i,2]] + gasPortOut.Ndot * x_v[mapping[i,1]] + reaction.Ndot[mapping[i,2]];
    end for;

//vdot_v=0.08;

  /* mole balance for inert substances (in reboiler only possible for liquids) */

  /* liquid phase */
    for i in 1:nSL loop
      if inert_Liquid[i] then
        der(V_liq*c_l[i]) = liquidPortIn.Ndot * x_l_in[i] + liquidPortOut.Ndot * x_l[i] + reaction.Ndot[i];
      end if;
    end for;

   HU_l = sum(c_l)*V_liq;
   HU_v = sum(c_v)*V_vap;

  vdot_l_in=F_in_l*MM_l_in/rho_l_in;
  vdot_l=F_out_l*MM_l/rho_l;
  vdot_v=F_out_v*MM_v/rho_v;

  /* energy balance */
  // der(HU_l*u_l+HU_v*u_v+m_col*cp_col*abs(T-T_ref))=F_in_l*h_l_in - F_out_l*h_l - F_out_v*h_v - Q_loss  + Q_in;
  der(HU_l*u_l+HU_v*u_v)=liquidPortIn.Ndot*h_l_in +liquidPortOut.Ndot*h_l +gasPortOut.Ndot*h_v - Q_loss  + Q_in + reaction.deltaH_R;

    /* summation equation */
    sum(x_l[i] for i in 1:nSL)=1;
    sum(x_v[i] for i in 1:nSV)=1;

    /* thermodynamic equilibrium */
    for i in 1:nS loop
    x_v[mapping[i,1]]= K[i] *x_l[mapping[i,2]];
    K[i] = thermoEquilibrium.K[i];
    end for;

checkMoleBal = F_in_l - F_out_l - F_out_v;

p_sat[:] = {mediumLiquid.p_sat[i] for i in 1:nSL};

// when time>20 then
//
// assert( vdot_v>=0, "Gasflow is negative - increase heat duty");
//
// end when;

algorithm

     for i in 1:nSL loop
        if inert_Liquid[i] then
           Ndot_l_transfer[i] := 0;
         else
           Ndot_l_transfer[i] :=-x_v[i]*F_out_v;
        end if;
     end for;

initial equation
if init_option==InitOptionEq.init_x then
  if init_standalone then

  T=T_init;
  eps_liq = eps_liq_init;

  for i in 1:nS-1 loop
  (HU_l*x_l[i]+HU_v*x_v[i])/(HU_l+HU_v)=x_total_start[i];
  end for;

  else

  T=T_init;
  p_sys = p_init;
  eps_liq = eps_liq_init;

  for i in 1:nS-1 loop
  (HU_l*x_l[i]+HU_v*x_v[i])/(HU_l+HU_v)=x_total_start[i];
  end for;

  end if;
elseif init_option==InitOptionEq.init_mol then
  if init_standalone then

  T=T_init;
  eps_liq = eps_liq_init;

  x_l[3]/(x_l[1]*0.018)=fixed_mol_init;

  else

  T=T_init;
  p_sys = p_init;
  eps_liq = eps_liq_init;

  x_l[3]/(x_l[1]*0.018)=fixed_mol_init;

  end if;
elseif init_option==InitOptionEq.init_heilbronn then
  if init_standalone then

  //T=T_init;
  // p_sys = p_init;
  //vdot_v=0.08;

  eps_liq =eps_liq_init;
  x_v[1]=x_total_start[1];
  x_l[3]/(x_l[1]*0.018)=fixed_mol_init;

  else

  //T=T_init;
  p_sys = p_init;
  //vdot_v=0.08;

  eps_liq =eps_liq_init;
  x_v[1]=x_total_start[1];
  x_l[3]/(x_l[1]*0.018)=fixed_mol_init;

  end if;
end if;
  annotation (DymolaStoredErrors,
    experiment,
    experimentSetupOutput,
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}),
         graphics),
    Diagram(graphics),
    Documentation(info="<html>
    <p><h4>Kettle reboiler - shell side</h4></p>
    <p>The reboiler is modelled as an evaporator with one liquid input and one liquid and vapour output respectively. The model consideres only the shell side of a kettle reboiler. On the inside the reboiler contains a weir which determines the maximum liquid level. Once the liquid level rises above the weir height a liquid stream starts to exit the reboiler. To model this process the francis-weir-formula is used.</p>
    <p>Mass and ernergy balances contain both liquid and vapour phases together. This means both phases are at thermodynamic equilibrum at any time.</p>
    <p>For now no startup processes from empty and cold states are considered. This means at t=0s the reboiler is already partially filled with liquid that has a temperature near the boiling point.</p>
    <p>So far no heat transfer resistance is considered for the heating. The amount of heat introduced into the reboiler is supplied by a <i>real Input</i> and directly accounted for in the energy balance. Also no pressure loss is calculated at the moment.</p>
    <p>The model can handle inert substances in the liquid phase. In case the liquid phase contains an inert the user has to specify which substances of the component vector are inert. If any chemical reactions occur they can also be taken into account. Reaction model and correlation for vapour-liquid equilibrium can both be selected via a dropdown menu.</p>
    <p>For initialisation the user can opt for either standalone operation or usage in larger systems because in standalone operation it is not possible to initialise the pressure since it is already given by the vapour sink.</p>
    </html>"));
end KettleReboilerEq;
