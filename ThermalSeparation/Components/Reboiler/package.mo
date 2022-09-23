within ThermalSeparation.Components;
package Reboiler "this package contains components which have both a liquid and a vapour phase in direct contact, with mass transfer between the two phases"
  extends Icons.Library.Red;

  model KettleReboilerEq_StartUpCCS

    import      Modelica.Units.SI;
    //import ThermalSeparation;
    extends ThermalSeparation.Icons.Color.Reboiler;
    outer ThermalSeparation.SystemTS systemTS;

  replaceable package MediumLiquid =
        ThermalSeparation.Media.H2O_CO2_MEA_Liq
    constrainedby ThermalSeparation.Media.BaseMediumLiquid
                                                 annotation(choicesAllMatching);
   replaceable model InnerHT =
      ThermalSeparation.HeatAndMassTransfer.HTResistance.NoHTResistance
      constrainedby ThermalSeparation.HeatAndMassTransfer.HTResistance.BaseHTResistance "heat transfer mechanism between bulk and wall"                     annotation(choicesAllMatching=true);

    InnerHT innerHT(n=1,T={T},A={A_HT},Qdot={Q_in},p={p_sys});

    parameter ThermalSeparation.Components.Reboiler.InitOptionEq
      init_option=ThermalSeparation.Components.Reboiler.InitOptionEq.init_x                                         annotation(Dialog(tab="Initialization"),Evaluate=true); // Enumerations.InitializationOption.init_x "initialization options"

    parameter Modelica.Units.SI.Area A_HT=5 "heat exchange area";

    replaceable package MediumVapour =
        ThermalSeparation.Media.H2O_CO2_Vap
    constrainedby ThermalSeparation.Media.BaseMediumVapour
                                                 annotation(choicesAllMatching);
    replaceable model ThermoEquilibrium =
        ThermalSeparation.PhaseEquilibrium.H2O_CO2_MEA_startUp
        constrainedby ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium annotation(choicesAllMatching);

    replaceable model Reaction =
        ThermalSeparation.Reaction.ReactionEquilibrium.CO2_MEA
    constrainedby ThermalSeparation.Reaction.BaseReaction
                                                annotation(choicesAllMatching);
    replaceable record Geometry =
         ThermalSeparation.Geometry.BasicGeometry annotation(choicesAllMatching);

  MediumVapour.BaseProperties mediumVapour(
        T0=T_ref, p=p_sys, c=c_vap, x=y, x_star=y,T=T);

  MediumLiquid.BaseProperties mediumLiquid(
        T0=T_ref, p=p_sys, x=x,T=T,h=h_liq);

  MediumLiquid.BaseProperties mediumLiquidIn(
        T0=T_ref, p=p_sys, x=x_in,T=T_liq_in,h=h_liq_in);

  ThermoEquilibrium thermoEquilibrium(
        redeclare replaceable package MediumVapour=MediumVapour,
        redeclare replaceable package MediumLiquid=MediumLiquid,
        nS=nS, mapping=mapping, T=T,x_l=x,x_v=y,v_v=MM_vap/rho_vap,p=p_sys,
        p_sat=p_sat,x_vap_liq=fill(1/nS,nS),startUp=startupCalc);

  Reaction reaction(
        redeclare replaceable package MediumLiquid=MediumLiquid,
        redeclare replaceable record Geometry=Geometry,
        propsLiq=
          mediumLiquid.properties,
        gamma=fill(1,nL),Ndot_l_transfer=N_liq_trans,c=c_liq, V=V_liq,n=1);

  //parameter

  //initialization
  parameter Real eps_liq_init = 0.2 annotation(Dialog(tab = "Initialization"));
  parameter SI.Temperature T_init = 300 annotation(Dialog(tab= "Initialization"));
  parameter SI.Pressure p_init = 2e5 annotation(Dialog(tab="Initialization"));
  parameter SI.MoleFraction x_init = 1/nS annotation(Dialog(tab="Initialization"));
  //total mole fractions substances
  parameter SI.MoleFraction x_total_start[nS]=fill(1/nS,nS) annotation(Dialog(tab="Initialization"));

  parameter Boolean initOption_TepsXfixed=true annotation(Dialog(tab="Initialization"));
  parameter Boolean initOption_xtotalfixed=false annotation(Dialog(tab="Initialization"));
  parameter Boolean initOption_withoutT=false annotation(Dialog(tab="Initialization"));
  parameter Boolean initOption_standalone=true annotation(Dialog(tab="Initialization"));
  parameter Boolean initOption_startup_xfixed=false annotation(Dialog(tab="Initialization"));
  parameter Boolean initOption_startup_xtotalfixed=false annotation(Dialog(tab="Initialization"));
  parameter Boolean initOption_startup_inert_Dyn=false annotation(Dialog(tab="Initialization"));
  parameter Boolean initOption_startup_RebplusDes=false annotation(Dialog(tab="Initialization"));
  parameter Boolean initOption_startup_RebplusDesalone=false annotation(Dialog(tab="Initialization"));
  parameter Boolean startup=false annotation(Dialog(tab="Start-up"));

  //inert in liquid phase
  parameter Boolean inert_liq[nL]={false,false,true};

  //number of substances which are equal in liquid and vapor phase
  parameter Integer nS=2;
  //number of substances in liquid phase
  final parameter Integer nL=MediumLiquid.nSubstance;
  //number of substances in vapor phase
  final parameter Integer nV=MediumVapour.nSubstance;

  //N_vap_out parameter k
  Real k=Vdot_nom/deltap_nom;
  parameter Real Vdot_nom=1;
  parameter Real deltap_nom=0.005;
  parameter Real beta_eps=0.001 annotation(Dialog(tab="Initialization"));

  //transfer over phase boundary
  SI.MolarFlowRate N_liq_trans[nL];

  //reaction
  SI.MolarFlowRate N_reaction[nL]=reaction.Ndot;
  Real delH_R=reaction.deltaH_R;

  //mapping
  parameter Integer mapping[nS,2]={{1,1},{2,2}};

  //Temperatures
  parameter SI.Temperature T_ref = systemTS.T_ref;
  SI.Temperature T(start=T_init);
  SI.Temperature T_liq_in(start=T_init);

  //pressure
  SI.Pressure p_sys;
  SI.Pressure p_sat[nL];
  SI.Pressure p_bub;

  //start-up
  parameter SI.Pressure p_amb=1e5;
  parameter SI.MoleFraction y_init[nS]=fill(1/nS,nS);
  parameter SI.Pressure p_friggel=0.3e5;
  Boolean startupCalc(start=startup);

  Real omega_eps;
  //Real omega_p;

  //Hold-ups
  SI.AmountOfSubstance HU_liq(start=100);
  SI.AmountOfSubstance HU_vap(start=100);

  //Vapor and liquid fraction
  Real eps_liq;
  Real eps_vap(min=0,max=1);

  //Spezifische innere Energie
  SI.MolarInternalEnergy u_liq=mediumLiquid.u;
  SI.MolarInternalEnergy u_vap=mediumVapour.u;

  //Volume Flow Rates
  SI.VolumeFlowRate V_liq_in;
  SI.VolumeFlowRate V_liq_out;
  SI.VolumeFlowRate V_vap_out(start=50);

  //Konzentrationen
  SI.Concentration c_liq_in[nL];
  SI.Concentration c_liq[nL];
  SI.Concentration c_vap[nV];

  Real dummy(stateSelect=StateSelect.prefer);
  Real dummy2(stateSelect=StateSelect.prefer);
  Real dummy3(stateSelect=StateSelect.prefer);

  //molar fractions of liquid and gas phase
  SI.MoleFraction y[nV](start=fill(1/nV,nV));
  SI.MoleFraction x[nL](start=fill(1/nL,nL));
  SI.MoleFraction x_in[nL];

  //density and molar mass of liquid phase
  SI.Density rho_liq=mediumLiquid.d;
  SI.MolarMass MM_liq=mediumLiquid.MM;

  SI.Density rho_liq_in=mediumLiquidIn.d;
  SI.MolarMass MM_liq_in=mediumLiquidIn.MM;

  //density and molar mass of vapor phase
  SI.Density rho_vap=mediumVapour.d;
  SI.MolarMass MM_vap(start=0.031)=mediumVapour.MM;

  //equilibrium constant
  Real K[nS](start=fill(1,nS));

  //Volume
  SI.Volume V_liq;
  SI.Volume V_vap;
  final SI.Volume V=A*H;
  SI.Height H_liq;

  //enthalpies
  SI.MolarEnthalpy h_liq_in;//=mediumLiquidIn.h;
  SI.MolarEnthalpy h_liq;//=mediumLiquid.h;
  SI.MolarEnthalpy h_vap=mediumVapour.h;

  //Heat flow rates
  SI.HeatFlowRate Q_in;
  parameter SI.HeatFlowRate Q_loss=0;
  SI.HeatFlowRate Q_loss_balance=if T>307 then Q_loss else 0;

  //MolarFlowRates

  SI.MolarFlowRate N_liq_in;
  SI.MolarFlowRate N_liq_out(start=1e-4);
  SI.MolarFlowRate N_vap_out(start=1e-4);
  SI.MolarFlowRate checkmolebalance;

  //geometry Reboiler
  parameter SI.Mass m_Reb=100;
  parameter SI.SpecificHeatCapacity c_p_Reb=0.46e3;
  parameter SI.Area A=1;
  parameter SI.Height H=1;
  parameter SI.Height H_weir=0.13;//0.6;
  parameter SI.Height L_weir=0.35;

    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation (
        Placement(transformation(extent={{76,-8},{96,12}}),
          iconTransformation(extent={{72,-10},{92,10}})));
    Interfaces.GasPortOut gasPortOut(redeclare package Medium = MediumVapour)
      annotation (Placement(transformation(extent={{-60,80},{-40,100}})));
    Interfaces.LiquidPortOut liquidPortOut(redeclare package Medium =
          MediumLiquid)
      annotation (Placement(transformation(extent={{-10,-100},{10,-80}})));
    Interfaces.LiquidPortIn liquidPortIn(redeclare package Medium =
        MediumLiquid)
      annotation (Placement(transformation(extent={{40,80},{60,100}})));
  equation

  dummy=c_liq[2];
  dummy2=HU_liq*u_liq + HU_vap*u_vap + m_Reb*c_p_Reb*abs(T-T_ref);
  dummy3=c_liq[1];

  //start-up

  //possibilities pressure-loss

  //1
  // if p_sys<gasPortOut.p then
  //   V_vap_out=0;
  // else
  //   p_sys=gasPortOut.p;
  // end if;

  //2
  //  if startupCalc then
  //    V_vap_out=max(0,k*(p_sys-gasPortOut.p));
  //  else
  //    p_sys=gasPortOut.p;
  //  end if;
  //p_sys=gasPortOut.p;
  //3
  V_vap_out=max(0,k*(p_sys-gasPortOut.p));
  //p_sys=gasPortOut.p;

  //ports
  gasPortOut.Ndot=-N_vap_out;
  gasPortOut.h_outflow=h_vap;
  gasPortOut.x_outflow=y;

  //V_vap_out=max(0,k*(p_sys-gasPortOut.p)); //pressure loss

  liquidPortOut.Ndot=-N_liq_out;
  liquidPortOut.h_outflow=h_liq;
  liquidPortOut.x_outflow=x;

  liquidPortIn.Ndot=N_liq_in;
  liquidPortIn.h_outflow=h_liq;
  inStream(liquidPortIn.h_outflow)=h_liq_in;
  inStream(liquidPortIn.x_outflow)=x_in;
  liquidPortIn.x_outflow=x;
  liquidPortIn.p=p_sys;

  heatPort.Q_flow = Q_in;
  heatPort.T = innerHT.Twall[1];
  der(HU_liq*u_liq + HU_vap*u_vap + m_Reb*c_p_Reb*abs(T-T_ref))= N_liq_in*h_liq_in-N_liq_out*h_liq-N_vap_out*h_vap + Q_in - Q_loss + delH_R;

  //energy balance

  //mole balance

  for i in 1:nS loop
     //der(V_vap*c_vap[i] + V_liq*c_liq[i]) = V_liq_in*c_liq_in[i] - V_vap_out*c_vap[i] - V_liq_out*c_liq[i];
    der(V_vap*c_vap[i] + V_liq*c_liq[mapping[i,2]]) = N_liq_in*x_in[mapping[i,2]] - N_vap_out*y[i] - N_liq_out*x[mapping[i,2]] + N_reaction[mapping[i,2]];
  end for;

  //inert

  for i in 1:nL loop
    if inert_liq[i] then
      der(V_liq*c_liq[i])=N_liq_in*x_in[i] - N_liq_out*x[i] + N_reaction[i];
    end if;
  end for;

  checkmolebalance=N_liq_in-N_liq_out-N_vap_out;

  //epsilon
  eps_vap=V_vap/V;
  eps_liq=V_liq/V;

   if startup then
     eps_vap=omega_eps*(1-eps_liq); //start up
     //eps_vap=1-eps_liq;
     omega_eps=smooth(2,if startupCalc then 1 + tanh(beta_eps*(p_bub-(p_sys-p_friggel))) else 1);
     //omega_eps=smooth(2,1 + tanh(beta_eps*(p_bub-(p_sys-p_friggel))));
   else
     eps_vap+eps_liq=1;
     omega_eps=0;
  end if;

  //equilibrium
  //     for i in 1:nS loop
  //       if startupCalc then
  //         y[i]=y_init[i];
  //       else
  //         y[mapping[i,1]]=K[i]*x[mapping[i,2]];
  //       end if;
  //         K[i]=thermoEquilibrium.K[i];
  //   end for;

  for i in 2:nS loop
       K[i]=thermoEquilibrium.K[i];
  end for;

  for i in 1:nS loop
    y[mapping[i,1]]= K[i]*x[mapping[i,2]];
  end for;

  //SaturationPressure
  for i in 1:nL loop
    p_sat[i]=mediumLiquid.p_sat[i];
  end for;

  //concentration --> molefraction
  for i in 1:nL loop
    c_liq[i]=x[i]*rho_liq/MM_liq;
  end for;
  for i in 1:nV loop
    c_vap[i]=y[i]*rho_vap/MM_vap;
  end for;
  for i in 1:nL loop
    c_liq_in[i]=x_in[i]*rho_liq_in/MM_liq_in;
  end for;

  //Hold-ups
  HU_liq=sum(c_liq)*V_liq;
  HU_vap=sum(c_vap)*V_vap;

  //Sums
  if startup then
    if startupCalc then
      p_sys=p_amb;
    else
      //y[mapping[1,1]]= K[1]*x[mapping[1,2]];
      K[1]=thermoEquilibrium.K[1];
     //sum(y[:])=1;
    end if;
  else
   //y[mapping[1,1]]= K[1]*x[mapping[1,2]];
    K[1]=thermoEquilibrium.K[1];
    //sum(y[:])=1;
  end if;

  sum(x[:])=1;
  sum(y[:])=1;
  //sum(x[i] for i in 1:nL)=1;

  //p_bub=sum(x.*p_sat);
  //p_bub=(x[1]/(x[1]+x[2]))*p_sat[1]+(x[2]/(x[1]+x[2]))*p_sat[2];
  p_bub=thermoEquilibrium.p_bubble;

  //height liquid
  H_liq=V_liq/A;

  //Francis-Weir-equation
  N_liq_out = noEvent(if H_liq<=H_weir then 0 else rho_liq/MM_liq*1.848*L_weir*(abs(H_liq-H_weir))^1.5);

  //mole --> volume
  N_liq_out=V_liq_out*rho_liq/MM_liq;
  N_liq_in=V_liq_in*rho_liq/MM_liq;
  N_vap_out=V_vap_out*rho_vap/MM_vap;

  // if p_bub<=p_sys then
  //    startupCalc=true;
  // else
  //    startupCalc=false;
  // end if;
  algorithm
  //transfer
  for i in 1:nL loop
    if inert_liq[i] then
      N_liq_trans[i]:=0;
    else
      N_liq_trans[i]:=-y[i]*N_vap_out;
    end if;
  end for;

  when omega_eps>=1 then
    startupCalc:=false;
    //reinit(y[1],0.9);
  end when;

  initial equation
    //N_liq_out=0;
    //p_sys=p_init;

    //eps_liq=0;
    //T=T_init;
    //y={0.5,0.5};
    //y[2]=0.5;
    //x[2]=x_init;
    //{{y[1]+y[2]=1;

    //T=T_init;
    //p_sys = p_init;
    //vdot_v=0.08;

    //eps_liq =eps_liq_init;
    //eps_vap=0;
    //y[1]=0.84;
    //x[3]/(x[1]*0.018)=7;

  if initOption_TepsXfixed then
    if initOption_standalone then
      eps_liq=eps_liq_init;
      T=T_init;
      if nL>2 or nV>2 then
       x[2]=x_init;
      end if;
    else
      eps_liq=eps_liq_init;
      p_sys=p_init;
      T=T_init;
      if nL>2 or nV>2 then
        x[2]=x_init;
      end if;
    end if;
  elseif initOption_xtotalfixed then
     if initOption_standalone then
       eps_liq=eps_liq_init;
       T=T_init;
       if nS==nL then
         for i in 1:nS-2 loop
           (HU_liq*x[i]+HU_vap*y[i])/(HU_liq+HU_vap)=x_total_start[i];
         end for;
       else
         for i in 1:nS-1 loop
           (HU_liq*x[i]+HU_vap*y[i])/(HU_liq+HU_vap)=x_total_start[i];
         end for;
       end if;
     else
       eps_liq=eps_liq_init;
       p_sys=p_init;
       T=T_init;
       if nS==nL then
         for i in 1:nS-2 loop
           (HU_liq*x[i]+HU_vap*y[i])/(HU_liq+HU_vap)=x_total_start[i];
         end for;
       else
         for i in 1:nS-1 loop
           (HU_liq*x[i]+HU_vap*y[i])/(HU_liq+HU_vap)=x_total_start[i];
         end for;
       end if;
     end if;
  elseif initOption_withoutT then
      if initOption_standalone then
       eps_liq=eps_liq_init;
       x[2]=x_init;
       if nS==nL then
         for i in 1:nS-2 loop
           (HU_liq*x[i]+HU_vap*y[i])/(HU_liq+HU_vap)=x_total_start[i];
         end for;
       else
         for i in 1:nS-1 loop
           (HU_liq*x[i]+HU_vap*y[i])/(HU_liq+HU_vap)=x_total_start[i];
         end for;
       end if;
      else
       eps_liq=eps_liq_init;
       p_sys=p_init;
       x[2]=x_init;
       if nS==nL then
         for i in 1:nS-2 loop
           (HU_liq*x[i]+HU_vap*y[i])/(HU_liq+HU_vap)=x_total_start[i];
         end for;
       else
         for i in 1:nS-1 loop
           (HU_liq*x[i]+HU_vap*y[i])/(HU_liq+HU_vap)=x_total_start[i];
         end for;
       end if;
      end if;
  elseif initOption_startup_xfixed then
      eps_liq=0;
      //eps_vap=0;
      T=T_init;
      x[2]=x_init;
  elseif initOption_startup_inert_Dyn then
      eps_liq=eps_liq_init;
      T=T_init;
      //p_sys=p_init;
  //       if inert_liq[2]==false then
  //         x[2]=x_init;
  //       else
  //         x[1]=x_init;
  //       end if;
      x[3]/(x[1]*0.018)=7;
      x[1]=x_init;
  elseif initOption_startup_xtotalfixed then
      eps_liq=eps_liq_init;
      //T=T_init;
        if nS==nL then
         for i in 1:nS-2 loop
           (HU_liq*x[i]+HU_vap*y[i])/(HU_liq+HU_vap)=x_total_start[i];
         end for;
       else
         for i in 1:nS-1 loop
           (HU_liq*x[i]+HU_vap*y[i])/(HU_liq+HU_vap)=x_total_start[i];
         end for;
        end if;
  elseif initOption_startup_RebplusDes or init_option==InitOptionEq.init_startUp then
    //y[2]=0.5;
    T=T_init;
    //Q_in=100e3;
    //p_sys = p_init;
    //eps_vap=0;
    eps_liq=eps_liq_init;
    x[3]/(x[1]*0.018)=7;
    x[1]=x_init;
  elseif initOption_startup_RebplusDesalone then
    //y[2]=0.5;
    T=T_init;
    //p_sys = p_init;
    //eps_vap=0;
    //eps_liq=eps_liq_init;
    x[3]/(x[1]*0.018)=7;
    x[1]=x_init;
  end if;
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}}),       graphics), Icon(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
  end KettleReboilerEq_StartUpCCS;

  type InitOptionEq = enumeration(
    init_x   "initializing using x_total",
    init_mol   "initializing using fixed molarity",
    init_heilbronn   "initialisation using p,eps,x_total[1] and molalityenumeration: initialization options for equilibrium reboiler",
    init_startUp   "initialization start-up",
    initOption_startup_inert_Dyn   "initialization start-up Dynstart") "enumeration: initialization options for equilibrium reboiler";
end Reboiler;
