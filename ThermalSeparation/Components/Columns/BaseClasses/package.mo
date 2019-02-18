within ThermalSeparation.Components.Columns;
package BaseClasses "base classes for columns"

  partial model BaseColumn_newStartUpShutDown
    "Gesamtmolbilanz, x, h, Ndot im Konnektor"

  outer ThermalSeparation.SystemTS systemTS;
  parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

    parameter Integer n(min=1)=1 "packed column: number of discrete elements in the section; plate column: number of trays in one section";

  //Initialization
    parameter SI.Pressure p_v_start_inlet = 1.9e5 annotation(Dialog(tab="Initialization"));
    parameter SI.Pressure p_v_start_outlet = 1.8e5 annotation(Dialog(tab="Initialization"));
    final parameter SI.Pressure p_v_start[n] = if n==1 then {p_v_start_outlet} else linspace(p_v_start_inlet,p_v_start_outlet,n);
    parameter Boolean x_l_profile = false  annotation(Dialog( tab="Initialization"));
    parameter Boolean x_v_profile = false annotation(Dialog( tab="Initialization"));
    parameter SI.MoleFraction x_l_start_const[nSL]= fill(1/nSL,nSL)                                 annotation(Dialog(enable=not x_l_profile, tab="Initialization"));
    parameter SI.MoleFraction x_v_start_const[nSV] = fill(1/nSV,nSV)      annotation(Dialog(enable=not x_v_profile, tab="Initialization"));
    parameter SI.MoleFraction x_l_start_in[nSL]= fill(1/nSL,nSL)                                 annotation(Dialog(enable=x_l_profile, tab="Initialization"));
    parameter SI.MoleFraction x_l_start_out[nSL]= fill(1/nSL,nSL)                                 annotation(Dialog(enable=x_l_profile, tab="Initialization"));
    parameter SI.MoleFraction x_v_start_in[nSV] = fill(1/nSV,nSV)         annotation(Dialog(enable=x_v_profile, tab="Initialization"));
    parameter SI.MoleFraction x_v_start_out[nSV] = fill(1/nSV,nSV)        annotation(Dialog(enable=x_v_profile, tab="Initialization"));
    final parameter SI.MoleFraction x_l_start[n,nSL](each fixed=false);
    final parameter SI.MoleFraction x_v_start[n,nSV](each fixed=false);
    parameter Real x_total_start[nSV]=fill(1/nSV,nSV) "total mole fraction in system (vapour and liquid), component ordering as in vapour medium"
                                                                                            annotation(Dialog(tab="Initialization"));

    parameter Boolean T_l_profile = false annotation(Dialog(tab="Initialization"));
    parameter Boolean T_v_profile = false annotation(Dialog(tab="Initialization"));
    parameter SI.Temperature T_vap_start_bottom = 300 annotation(Dialog(tab="Initialization"));
    parameter SI.Temperature T_vap_start_top = 300 annotation(Dialog(tab="Initialization"));
    parameter SI.Temperature T_liq_start_bottom = 300 annotation(Dialog(tab="Initialization"));
    parameter SI.Temperature T_liq_start_top = 300 annotation(Dialog(tab="Initialization"));
    parameter SI.Temperature T_vapour_start= 300 annotation(Dialog(tab="Initialization"));
    parameter SI.Temperature T_liquid_start= 300 annotation(Dialog(tab="Initialization"));
   //for equilibrium models, vapour start temperature always equals liquid start temperature
  //   final parameter SI.Temperature T_v_start[n]= if EQ then (if T_l_profile then linspace(T_liq_start_bottom, T_liq_start_top, n) else ones(n)*T_liquid_start) else (if T_v_profile and not n==1 then linspace(T_vap_start_bottom, T_vap_start_top, n) else (if T_v_profile and n==1 then ones(n)*(T_vap_start_bottom+ T_vap_start_top)/2 else ones(n)*T_vapour_start));
  //   final parameter SI.Temperature T_l_start[n]= if EQ then (if T_l_profile then linspace(T_liq_start_bottom, T_liq_start_top, n) else ones(n)*T_liquid_start) else (if T_l_profile and not n==1 then linspace(T_liq_start_bottom, T_liq_start_top, n) else (if T_l_profile and n==1 then ones(n)*(T_liq_start_bottom+ T_liq_start_top)/2 else ones(n)*T_liquid_start));
    final parameter SI.Temperature T_v_start[n] = if (T_v_profile and not n==1) then linspace(T_vap_start_bottom, T_vap_start_top, n) else (if (T_v_profile and n==1) then ones(n)*(T_vap_start_bottom+ T_vap_start_top)/2 else ones(n)*T_vapour_start);
    final parameter SI.Temperature T_l_start[n] = if (T_l_profile and not n==1) then linspace(T_liq_start_bottom, T_liq_start_top, n) else (if (T_l_profile and n==1) then ones(n)*(T_liq_start_bottom+ T_liq_start_top)/2 else ones(n)*T_liquid_start);

    // Homotopy-Method
    replaceable model HomotopyMethod =
      ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.NoHomotopy
                                  constrainedby
    ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.BaseHomotopy
                                                                                                                                      "activate homotopy and set nominal values"
                                                                            annotation(Dialog(tab="Initialization", group="Homotopy"),choicesAllMatching=true);

  HomotopyMethod homotopyMethod(nS=nS,n=n,nSV=nSV,nSL=nSL);

  //Result-Record
  Results results(n=n, nSL=nSL, nSV=nSV, T_l=T_l, x_l=x_l, x_v=x_v,  x_l_star=x_l_star, x_v_star=x_v_star,c_l = c_l, c_v=c_v, p_v=p_v, Vdot_l=Vdot_l, Vdot_v=Vdot_v, eps_liq = eps_liq, startUp=startUp);

  //Equilibrium or Non-Equilibrium model?
    parameter Boolean EQ= false "equilibrium model is used, no mass transfer, value provided by film model"
                                                                                  annotation(Dialog(enable=false));
  //Medium
  parameter Integer mapping[nS,2] = {{i,i} for i in 1:nS} "parameter to map the different medium vectors one to another";
  parameter Boolean inertVapour[nSV] = fill(false,nSV) "true for each component which is inert in the vapour phase";
  parameter Boolean inertLiquid[nSL] = fill(false,nSL) "true for each component which is inert in the liquid phase";
  replaceable package MediumVapour = ThermalSeparation.Media.H2O_CO2_Vap
      constrainedby ThermalSeparation.Media.BaseMediumVapour "medium to be used in vapour phase"                                                         annotation(choicesAllMatching);
   MediumVapour.BaseProperties mediumVapour[n](c=c_v, each T0=T_ref,   p=p_v[1:n], T=T_v, x=x_v,  x_star=x_v_star);
   MediumVapour.BaseProperties mediumVapourIn(c=c_v_in, T0=T_ref, p=upStreamIn.p, T=T_v_in, x=x_v_in, x_star=x_v_in);
    MediumVapour.BaseProperties mediumVapourStar[n](c=c_v, each T0=T_ref, p=p_v[1:n], T=T_star, x=x_v_star, x_star=x_v_star);

  replaceable package MediumLiquid = ThermalSeparation.Media.H2O_CO2_MEA_Liq
     constrainedby ThermalSeparation.Media.BaseMediumLiquid "medium to be used in liquid phase"                                                         annotation(choicesAllMatching);
         MediumLiquid.BaseProperties mediumLiquid[n](each T0=T_ref,  p=p_hyd[1:n], T=T_l, x= x_l, h=if homotopyMethod.bool_h and homotopyMethod.useHomotopy then homotopy(actual=h_l,simplified=fill(homotopyMethod.h_liq,n)) else h_l);
   MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref,  p=p_hyd[n+1], T=T_l_in, x=x_l_in, h=h_l_in);
   MediumLiquid.BaseProperties mediumLiquidStar[n](each T0=T_ref, p=p_hyd[1:n], T=T_star, x=x_l_star,h=h_l_star);
   MediumLiquid.ActivityCoefficient activityCoeff[n](T=T_star,x_l=x_l_star);
   MediumVapour.EvaporationEnthalpy evapEnthalpy[n](  p=p_hyd[1:n], T=T_v);
   constant Boolean h_evap_medium = MediumVapour.delta_hv_medium;
   ThermalSeparation.Units.MolarEnthalpy delta_hv[n,nSV] = if h_evap_medium then zeros(n,nSV) else evapEnthalpy.h;

    parameter Integer nS(min=2) "number of species which are equal in vapour and liquid phase";
    final parameter Integer nL=MediumLiquid.nSubstance -nS "number of additional substances which are only in liquid phase";
    final parameter Integer nV = MediumVapour.nSubstance -nS "number of additional substances which are only in the vapour phase";
    final parameter Integer nSL = MediumLiquid.nSubstance;
    final parameter Integer nSV = MediumVapour.nSubstance;

    /*** vapour properties ***/
    SI.Density rho_v[n]= if homotopyMethod.bool_rho and homotopyMethod.useHomotopy then homotopy(actual=mediumVapour.d,simplified=fill(homotopyMethod.rho_vap,n)) else mediumVapour.d "mixture vapour density";
    // SI.Density rho_v[n]= mediumVapour.d "mixture vapour density";
    SI.Density rho_v_in =  mediumVapourIn.d;
    SI.MolarMass MM_v[n]( start=0.028*ones(n), each stateSelect=StateSelect.prefer)= mediumVapour.MM "molar mass of the vapour mixture ";
    SI.MolarMass MM_v_in( start=0.03) = mediumVapourIn.MM;
    ThermalSeparation.Units.MolarEnthalpy h_v[n] = if homotopyMethod.bool_h and homotopyMethod.useHomotopy then homotopy(actual=mediumVapour.h,simplified=fill(homotopyMethod.h_vap,n)) else mediumVapour.h;
    //ThermalSeparation.Units.MolarEnthalpy h_v[n] = mediumVapour.h;
    ThermalSeparation.Units.MolarEnthalpy h_v_in = mediumVapourIn.h;
    SI.MolarInternalEnergy u_v[n](each stateSelect=StateSelect.default)= mediumVapour.u;
    SI.Concentration c_v_star[n,nSV];
    SI.Density rho_v_star[n] = mediumVapourStar.d;

    /*** liquid properties ***/
    SI.Density rho_l[n]= if homotopyMethod.bool_rho and homotopyMethod.useHomotopy then homotopy(actual=mediumLiquid.d,simplified=fill(homotopyMethod.rho_liq,n)) else mediumLiquid.d "mixture liquid density";
    //   SI.Density rho_l[n] = mediumLiquid.d "mixture liquid density";
    SI.Density rho_l_in = mediumLiquidIn.d;
    SI.Density rho_l_star[n]= mediumLiquidStar.d;
    SI.MolarMass MM_l[n](start=fill(0.018,n))= mediumLiquid.MM "molar mass of the liquid mixture";
    SI.MolarMass MM_l_in= mediumLiquidIn.MM;
    SI.MolarMass MM_l_star[n]= mediumLiquidStar.MM;
    ThermalSeparation.Units.MolarEnthalpy h_l[n];//=if homotopyMethod.bool_h and homotopyMethod.useHomotopy then homotopy(actual=mediumLiquid.h,simplified=fill(homotopyMethod.h_liq,n)) else mediumLiquid.h;
    //ThermalSeparation.Units.MolarEnthalpy h_l[n];
    ThermalSeparation.Units.MolarEnthalpy h_l_in;//=mediumLiquidIn.h;
    ThermalSeparation.Units.MolarEnthalpy h_l_star[n];
    SI.MolarInternalEnergy u_l[n](each stateSelect=StateSelect.default) =  mediumLiquid.u;

  //Variables upStream
    SI.Concentration c_v_in[nSV];
    SI.Concentration c_v[n,nSV](each stateSelect=StateSelect.default) annotation(Dialog(group="Initialization",showStartAttribute=true));
    SI.MoleFraction x_v_in[nSV];
    parameter SI.MoleFraction x_v_dummy[nSV]={1,0};
    SI.MoleFraction x_v[n,nSV](start=x_v_start);
    SI.VolumeFlowRate Vdot_v_in(nominal=1e-2);
    SI.VolumeFlowRate Vdot_v[n](nominal=fill(1e-2,n));
    SI.Temperature T_v_in;

    SI.Pressure p_v[n+1]( start=fill(1.0e5,n+1)) "p_v[j] = pressure on the j-th stage, p_v[n+1] is the pressure in the first element of the sucesseding component";
    SI.Temperature T_v[n](nominal=fill(350,n),start=fill(350,n));

    //Variables downStream
    SI.Concentration c_l_in[nSL] "molar concentration in the liquid at the liquid outlet of each stage";
    SI.Concentration c_l[n,nSL](each stateSelect=StateSelect.default) annotation(Dialog(group="Initialization",showStartAttribute=true));

    SI.MoleFraction x_l_in[nSL];
    SI.MoleFraction x_l[n,nSL](start=x_l_start);
    SI.VolumeFlowRate Vdot_l_in(start=0.03, nominal=1e-4);
    SI.VolumeFlowRate Vdot_l[n](nominal=fill(1e-4,n));
    SI.Temperature T_l_in;
    SI.Temperature T_l[n];
    SI.Concentration c_l_star[n,nSL];

    //Model for reaction: to be supplied by extending class
    SI.MolarFlowRate Ndot_reac[n,nSL];
    SI.HeatFlowRate Qdot_reac[n];

    //instances of connectores to next stage
    ThermalSeparation.Interfaces.GasPortIn upStreamIn(redeclare package Medium =
          MediumVapour) annotation (Placement(transformation(extent={{-80,-100},{-60,-80}},
                         rotation=0), iconTransformation(extent={{-80,-100},{-60,-80}})));
    ThermalSeparation.Interfaces.GasPortOut upStreamOut(redeclare package Medium =
          MediumVapour) annotation (Placement(transformation(extent={{-80,80},{-60,
              100}}, rotation=0), iconTransformation(extent={{-80,80},{-60,100}})));
    ThermalSeparation.Interfaces.LiquidPortIn downStreamIn(redeclare package
        Medium = MediumLiquid) annotation (Placement(transformation(extent={{60,
              80},{80,100}}, rotation=0), iconTransformation(extent={{60,80},{80,
              100}})));
    ThermalSeparation.Interfaces.LiquidPortOut downStreamOut(redeclare package
        Medium = MediumLiquid) annotation (Placement(transformation(extent={{60,-100},
              {80,-80}}, rotation=0), iconTransformation(extent={{60,-100},{80,-80}})));

    //initial equation for eps_liq is supplied in the extending class!
    SI.VolumeFraction eps_liq[n](each stateSelect=StateSelect.default) "liquid volume fraction";
    SI.VolumeFraction eps_vap[n](start=fill(0.99,n)) "vapour volume fraction";
    SI.Temperature T[n];
    SI.HeatFlowRate Qdot_wall[n] "heat flow rate to wall";

    //to be supplied in extending class
    SI.MolarFlowRate Ndot_v_transfer[        n,nSV](start=fill(-0.1,n,nSV));
    SI.MolarFlowRate Ndot_l_transfer[        n,nSL](start=fill(0.1,n,nSL));

      SI.HeatFlowRate Edot_l_transfer[n];
    SI.HeatFlowRate Edot_v_transfer[n];

    SI.Temperature T_star[n](start=fill(293.15,n));

  SI.Pressure p_v_in(stateSelect=StateSelect.default,start=1.5e5);

  SI.Pressure p_sat[n,nSL];
  SI.Pressure p_sat_bulk[n,nSL];

  /*** Feed variables: has to be supplied by the extending class ***/
  SI.VolumeFlowRate Vdot_v_feed[n];
  SI.Concentration c_v_feed[n,nSV];
  SI.SpecificEnthalpy h_v_feed[n];
  SI.VolumeFlowRate Vdot_l_feed[n];
  SI.Concentration c_l_feed[n,nSL];
  SI.SpecificEnthalpy h_l_feed[n];
  SI.Density rho_l_feed[n];
  SI.Density rho_v_feed[n];
  SI.MolarMass MM_l_feed[n];
  SI.MolarMass MM_v_feed[n];

     SI.MoleFraction x_l_star[n,nSL](start=x_l_start);
   SI.MoleFraction x_v_star[n,nSV](start=x_v_start);

  /*** thermodynamic equilibrium ***/
   replaceable model ThermoEquilibrium =
       ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid                                  constrainedby
    ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium                                                                                                                  "model for phase equilibrium"
                                                                                            annotation(choicesAllMatching=true);
      ThermoEquilibrium bubblePressure[n](each nS=nS,  each mapping =                            mapping,
      redeclare replaceable package MediumVapour =   MediumVapour,
    redeclare replaceable package MediumLiquid =
       MediumLiquid, p=ones(n)*p_initial, T=T_l, x_v=x_v, x_l=x_l, p_sat=p_sat_bulk, v_v=MM_v./rho_v,each x_vap_liq=fill(1,nS));

protected
    parameter SI.Area A "cross sectional area of the column";
    parameter SI.Length H "height of this section of the column";
    parameter Real eps "void fraction in the column";
    final parameter SI.Mass mass_solid[n]=A*H/n*(1 - eps)*
        rho_solid;
    parameter SI.Density rho_solid[n];
    parameter SI.SpecificHeatCapacity c_solid;

   /***variables for reorganisation of mediums in the medium vector ***/
  //parameter Integer mL( fixed=false);
  // parameter Integer mV(fixed=false);

public
    Boolean bool_eps[n];

  /*** entrainment ***/
    SI.VolumeFlowRate Vdot_le[n] "liquid volume flow entrained by vapour"; // has to be supplied in extending class

  /*** StartUp ***/
    parameter Boolean considerStartUp = false "true if StartUp is to be considered"
                                            annotation(Dialog(tab="StartUp"));
    parameter Boolean considerShutDown=false "true if ShutDown is to be considered"  annotation(Dialog(tab="ShutDown"));
    parameter Boolean StartUp_CCS=false "true if StartUp of carbon capture plant is to be considered"
                                                        annotation(Dialog(tab="StartUp",enable=considerStartUp));
    parameter Boolean switchingCondition_Boiling=true "true if boiling state is switching condition" annotation(Dialog(tab="StartUp",enable=StartUp_CCS));
    parameter Boolean switchingCondition_Absorber_x_v=false "true if vapour composition is switching condition"
                                                        annotation(Dialog(tab="StartUp"),enable=StartUp_CCS);

    parameter Real x_v_switch=0.05 "vapour mole fraction value which is to be achieved"
                                                        annotation(Dialog(tab="StartUp",enable=switchingCondition_Absorber_x_v));
    parameter Integer componentNumber=3 "number of vapour component number in model"
                                                 annotation(Dialog(tab="StartUp",enable=switchingCondition_Absorber_x_v));
    parameter SI.Pressure p_initial= 1e5 "pressure in column at start-up" annotation(Dialog(tab="StartUp"));//ujpStreamOut.p ;
    parameter Real gain= 0.01 "controler gain to maintain initial pressure before switch" annotation(Dialog(tab="StartUp"));
    parameter Real k=0.2e-3 "large value for steep omega" annotation(Dialog(tab="StartUp"));
    parameter Real friggelfaktor = 0.0002e5 "value added to p_bub to shift switching time" annotation(Dialog(tab="StartUp"));
    parameter Boolean smooth_startUp=false "true if smooth switching is to be considered"
                                                     annotation(Dialog(tab="StartUp",group="Smooth Start-Up"));
    parameter Real delay_startUp= 200 "time delay for smooth startUp"
                                          annotation(Dialog(tab="StartUp",group="Smooth Start-Up",enable=smooth_startUp));

    parameter Boolean lowBoilingPoint[nSV]=fill(false,nSV) "true if substance has low boiling point"
                                                                                annotation(Dialog(tab="StartUp"));
    parameter Real y_PID=10 "maximal value for supply startUp PID controller" annotation(Dialog(tab="StartUp"));
    parameter Real Vdot_startUp_pressure=0.005 "value when supply PID controller is switched off" annotation(Dialog(tab="StartUp"));

    SI.Pressure p_bub[n]= bubblePressure.p_bubble "mixture bubble pressure";
    SI.Pressure p_hyd[n+1] "hydraulic pressure";
    Real omega[n];
    Boolean startUp[n](start=fill(true,n), each fixed=false);
    Boolean before_transition[n](start=fill(true,n),each fixed=false);
    Real transition_time[n](start=fill(1e7,n));
    //Real switch_time[n](start=fill(1e7,n));
    Real Ndot_source_startUp[n] "dummy molar flow rate to account for discharge of inert gas during startUp";
    //Real Ndot_source_startUp_l[n];
    Real Ndot_source_shutDown[n] "dummy molar flow rate to avoid stiff system";
    Real N_dummyShutDown=0.001 "dummy molar flow rate to avoid stiff system";
    //Real Ndot_source_shutDown_l[n] "dummy molar flow rate to avoid stiff system";
     //Real N_dummyShutDown_l=1 "dummy molar flow rate to avoid stiff system";

  /*** monitoring ***/
    //only for monitoring purpose: is the sum of the x really equal to one?
   Real sum_xl[n] = sum(x_l[:,i] for i in 1:nSL);
   Real sum_xv[ n] = sum(x_v[:,i] for i in 1:nSV);
   Real sum_xl_star[ n] = sum(x_l_star[:,i] for i in 1:nSL);
   Real sum_xv_star[ n] = sum(x_v_star[:,i] for i in 1:nSV);
   SI.MolarFlowRate Ndot_trans[nSL] = sum(Ndot_l_transfer[j,:] for j in 1:n);
     SI.MolarFlowRate Ndot_trans_vap[nSV] = sum(Ndot_v_transfer[j,:] for j in 1:n);

  Real Edot_l=sum(Edot_l_transfer);
  Real Edot_v=sum(Edot_v_transfer);
  SI.MassFlowRate mdot_v[n] = Vdot_v.*rho_v;
  SI.MassFlowRate mdot_l[n] = Vdot_l.*rho_l;
  Real X_v[n,nSV] "mass fraction vapour";
  Real X_v_in[nSV] "mass fraction vapour in";
  Real X_l[n,nSL] "mass fraction liquid";
  SI.Volume V_liq = sum(A*H/n*eps*eps_liq);
  SI.MolarFlowRate Ndot_v[n] "total molar flow rate vapour";
  SI.MolarFlowRate Ndot_v_in "total molar flow rate vapour";
  SI.MolarFlowRate Ndot_l[n] "total molar flow rate liquid";
  SI.MolarFlowRate Ndot_l_in "total molar flow rate vapour";
         Real n_i_liq[n,nSL];
     Real n_i_vap[n,nSV];
     Real n_liq[n] = sum(n_i_liq[:,i] for i in 1:nSL);
       Real n_vap[n] = sum(n_i_vap[:,i] for i in 1:nSV);
       Real n_total[n] = n_liq + n_vap;
  Real n_mol_L[n];
  Real n_mol_V[n];

  SI.MolarMass MM_v_star[n](start=fill(0.02,n))=mediumVapourStar.MM;
  Real check[nS-1];

  // ThermalSeparation.Utilities.LimPID_Input PID[
  //                                           n](each u_s = p_initial,u_m = if considerStartUp then p_v[1:n] else fill(p_initial,n),
  // yMax=100,
  // k=0.005,
  // initType=Modelica.Blocks.Types.InitPID.InitialOutput)   annotation (Placement(transformation(extent={{-24,-12},{-4,8}})));

          ThermalSeparation.Utilities.LimPID_Input PID[
                                                n](each u_s = p_initial,
      each controllerType=Modelica.Blocks.Types.SimpleController.P,
      each k=gain,
      each initType=Modelica.Blocks.Types.InitPID.NoInit,
      each yMax=y_PID)                                              annotation (Placement(transformation(extent={{-20,20},
              {0,40}})));

    Modelica.Blocks.Interfaces.BooleanInput StartUp_signal
      annotation (Placement(transformation(extent={{-102,-66},{-92,-56}})));
    Modelica.Blocks.Interfaces.BooleanInput ShutDown_signal
      annotation (Placement(transformation(extent={{-102,-50},{-92,-40}})));
    ThermalSeparation.Components.Columns.BaseClasses.InternalStepSequence_StateGraph1 internalStepSequence[n]
      annotation (Placement(transformation(extent={{-20,0},{0,20}})));
    Modelica.Blocks.Routing.BooleanReplicator shutDownReplicator(nout=n)
      annotation (Placement(transformation(extent={{-24,-48},{-16,-42}})));
    Modelica.Blocks.Routing.BooleanReplicator startUpReplicator(nout=n)
      annotation (Placement(transformation(extent={{-72,-64},{-62,-58}})));
    Modelica.Blocks.Sources.BooleanExpression booleanExpression(y=Vdot_v_in >
          Vdot_startUp_pressure)
      annotation (Placement(transformation(extent={{-64,6},{-50,18}})));
    Modelica.Blocks.Routing.BooleanReplicator endPhase1(nout=n)
      annotation (Placement(transformation(extent={{-38,8},{-30,14}})));
  initial algorithm

  /*** linear concentration profile at t=0 ***/
  if n==1 then
    for i in 1:nSV loop
      x_v_start[:,i] := if x_v_profile then fill(x_v_start_in[i],n) else fill(x_v_start_const[i],n);
    end for;
    for i in 1:nSL loop
      x_l_start[:,i] := if x_l_profile then fill(x_l_start_in[i],n) else fill(x_l_start_const[i],n);
    end for;
  else
  for j in 1:n loop
    for i in 1:nSV loop
      if x_v_profile then
        x_v_start[j,i] :=(x_v_start_out[i] - x_v_start_in[i])/(n - 1)*j +
            x_v_start_in[i] - (x_v_start_out[i] - x_v_start_in[i])/(n - 1);
      else
        x_v_start[j,i] :=x_v_start_const[i];
      end if;
      end for;
        for i in 1:nSL loop
      if x_l_profile then
        x_l_start[j,i] :=(x_l_start_out[i] - x_l_start_in[i])/(n - 1)*j +
            x_l_start_in[i] - (x_l_start_out[i] - x_l_start_in[i])/(n - 1);
      else
        x_l_start[j,i] :=x_l_start_const[i];
      end if;
      end for;
      end for;
  end if;

  equation
        for j in 1:n loop
       // PID[j].u_m = if considerStartUp  then p_v[j]  else p_initial;
              PID[j].u_m = if (considerStartUp and startUp[j]) then p_v[j] else  p_initial;

        end for;
    /***monitoring ***/
    for j in 1:n loop
    X_v[j,:]=x_v[j,:]/MM_v[j].*MediumVapour.MMX[:];
    X_l[j,:]=x_l[j,:]/MM_l[j].*MediumLiquid.MMX[:];
      n_i_liq[j,:]=c_l[j,:].*eps_liq[j]*A*H/n*eps;
      n_i_vap[j,:]=c_v[j,:].*eps_vap[j]*A*H/n*eps;
    end for;
    X_v_in[:] = x_v_in[:]/MM_v_in .* MediumVapour.MMX[:];
    //Verknpfung der Konnektorgren mit Variablen
    //upstream
    x_v_in = inStream(upStreamIn.x_outflow);
    h_v_in = inStream(upStreamIn.h_outflow);
    c_v_in = x_v_in / MM_v_in *rho_v_in;

    x_l_in = inStream(downStreamIn.x_outflow);
    h_l_in = inStream(downStreamIn.h_outflow);
    c_l_in = x_l_in / MM_l_in *rho_l_in;

    upStreamOut.x_outflow = x_v[n,:];
    upStreamOut.h_outflow = h_v[n];
    upStreamIn.x_outflow = x_v[1,:];
    upStreamIn.h_outflow = h_v[1];

    downStreamOut.x_outflow = x_l[1,:];
    downStreamOut.h_outflow = h_l[1];
    downStreamIn.x_outflow = x_l[n,:];
    downStreamIn.h_outflow = h_l[n];

    upStreamIn.p = p_v_in;
    upStreamOut.p = p_hyd[n+1];
      if ShutDown_signal and time >=100 then
        if nSV>=3 then
        N_dummyShutDown = Ndot_v_in;
         //upStreamIn.Ndot + N_dummyShutDown = Ndot_v_in;
        else
        upStreamIn.Ndot + N_dummyShutDown = Ndot_v_in;
        end if;
      else
       upStreamIn.Ndot = Ndot_v_in;
      end if;

  //     if ShutDown_signal and time >100 then
  //       if nSV>=3 then
  //       //downStreamIn.Ndot+ N_dummyShutDown_l = Ndot_l_in;
  //       N_dummyShutDown_l = Ndot_l_in;
  //       else
  //       downStreamIn.Ndot + N_dummyShutDown_l = Ndot_l_in;
  //       end if;
  //     else
  //      downStreamIn.Ndot = Ndot_l_in;
  //     end if;
    //upStreamIn.Ndot = Ndot_v_in;

     upStreamOut.Ndot = -Ndot_v[n];
   //upStreamOut.p_medium = p_hyd[n];
    downStreamIn.p = p_hyd[n+1];
    downStreamIn.Ndot = Ndot_l_in;
    downStreamOut.Ndot = -Ndot_l[1];

  Vdot_v= Ndot_v ./ rho_v .*MM_v "total molar flow rate vapour";
  Vdot_v_in= Ndot_v_in / rho_v_in*MM_v_in "total molar flow rate vapour";
  Vdot_l= Ndot_l ./ rho_l .*MM_l "total molar flow rate liquid";
  Vdot_l_in= Ndot_l_in / rho_l_in*MM_l_in "total molar flow rate vapour";

    /***eps_liq > 1: Betriebsendpunkt erreicht, allerdings ist beim Anfahrvorgang u.U. eps_liq = 1 ***/
    for j in 1:n loop
   // assert(eps_liq[j] < 0.99, "liquid volume fraction eps_liq gets larger than 1, the column can not be operated at this point!");
    end for;


  /*** StartUp ***/
  if considerStartUp then
    if StartUp_CCS then
      if switchingCondition_Boiling then
        //Desorber
        for j in 1:n loop
        if internalStepSequence[j].ColdState.active then
          p_hyd[j] = p_v[j];
          omega[j] = 1;
          Ndot_source_startUp[j]=PID[j].y;
          Ndot_source_shutDown[j]=0;
        elseif internalStepSequence[j].HeatUpState.active then
          p_hyd[j] = p_v[j];
          omega[j] = 1;
          Ndot_source_startUp[j]=0;
          Ndot_source_shutDown[j]=0;
        elseif internalStepSequence[j].ShutDown.active then
          p_hyd[j] = p_v[j];
          omega[j] = 1;
          if j==1 then
             Ndot_source_startUp[j]=0.0001;
             Ndot_source_shutDown[j]= - Ndot_source_startUp[j] - N_dummyShutDown;
          else
             Ndot_source_startUp[j]=0.0001;
             Ndot_source_shutDown[j]=-Ndot_source_startUp[j];
          end if;
        else
          p_hyd[j] = p_v[j];
          omega[j] = 1;
          Ndot_source_startUp[j]=0;
          Ndot_source_shutDown[j]=0;
        end if;
        end for;
      else
             //Absorber
        for j in 1:n loop
        if internalStepSequence[j].ColdState.active then
          p_hyd[j] = p_v[j];
          omega[j] = 1;
          Ndot_source_startUp[j]=PID[j].y;
          Ndot_source_shutDown[j]=0;
        elseif internalStepSequence[j].ShutDown.active then
          p_hyd[j] = p_v[j];
          omega[j] = 1;
          if j==1 then
             Ndot_source_startUp[j]=0.0001;
             Ndot_source_shutDown[j]= - Ndot_source_startUp[j] - N_dummyShutDown;
          else
             Ndot_source_startUp[j]=0.0001;
             Ndot_source_shutDown[j]=-Ndot_source_startUp[j];
          end if;
        else
          p_hyd[j] = p_v[j];
          omega[j] = 1;
          Ndot_source_startUp[j]=0;
          Ndot_source_shutDown[j]=0;
        end if;
        end for;
      end if;
     else
         for j in 1:n loop
            when p_bub[j] + friggelfaktor >= p_initial then // der Druck, den auch das Medienmodell sieht, friggelfaktor, weil sonst Umschalten zu spt
              startUp[j] = false;
            end when;
  //      if p_bub[j] <= p_v[j] then
  //        omega[j] = min(1,1 + tanh(k*(p_bub[j] + friggelfaktor - p_initial)));
  //      else
  //        omega[j] = 1;
  //      end if;
            omega[j] =  min(1,1 + tanh(k*(p_bub[j]+ friggelfaktor   - p_initial)));
   //   p_hyd[j] = if startUp[j] then p_initial else p_v[j];
            p_hyd[j] =  (p_initial * (1 - omega[j])) + p_v[j] * omega[j];

            if startUp[j] then
             Ndot_source_startUp[j]=PID[j].y;
            else
             Ndot_source_startUp[j]=0;
            end if;
          end for;
    end if;
  else  startUp = fill(false,n);
        for j in 1:n loop Ndot_source_startUp[j]=0; end for;
        omega = ones(n);
        p_hyd = p_v;
        before_transition = fill(false,n);
        transition_time=fill(0,n);
        startUp = fill(false,n);
  end if;

  p_hyd[n+1] = p_v[n+1];

  /*** give startup information to film model ***/
  for j in 1:n loop
    if internalStepSequence[j].ColdState.active or internalStepSequence[j].HeatUpState.active then
      before_transition[j]=true;
    else
      before_transition[j]=false;
    end if;
    when before_transition[j]==false then
      transition_time[j] = time;
    end when;
    when time>=(transition_time[j]+2*delay_startUp) and internalStepSequence[j].NormalOperation.active then
            startUp[j] = false;
    end when;
  end for;

  /*** calculation of molar fraction using the pressure and the molar concentration ***/
  for j in 1:n loop

      c_l_star[j,:] =x_l_star[j,:] / MM_l_star[j]*rho_l_star[j];
      c_v_star[j,:] =x_v_star[j,:] / MM_v_star[j]*rho_v_star[j];
      // c_l_star[j,i] = (x_l_star[j,i]+x_l[j,i])/2 / MM_l_star[j]*rho_l_star[j];
    for i in 1:nS loop
      // x_l[j,i] = if use_v then c_l[j,i] * mediumLiquid[j].v else c_l[j,i] * MM_l[j] /rho_l[j];
      x_l[j,i] = c_l[j,i]  * mediumLiquid[j].v;
      x_v[j,i] = c_v[j,i] * MM_v[j] /rho_v[j];
    end for;
    for i in 1:nL loop
      x_l[j,i+nS] =  c_l[j,i+nS] *mediumLiquid[j].v;
      //x_l[j,i+nS] = if use_v then c_l[j,i+nS] *mediumLiquid[j].v else c_l[j,i+nS] * MM_l[j] /rho_l[j];
    end for;
    for i in 1:nV loop
      x_v[j,i+nS] = c_v[j,i+nS] * MM_v[j] /rho_v[j];
    end for;
    end for;
    for i in 1:nSL loop

    end for;

    if n == 1 then
       bool_eps[1] = if (eps_liq[1]<1e-5 and Vdot_l_in<1e-8) then true else false;
     else
       bool_eps[1] = if eps_liq[1]<1e-5 and Vdot_l[2]<1e-8 then true else false;

    for j in 2:n-1 loop
       bool_eps[j] = if eps_liq[j]<1e-5 and Vdot_l[j+1]<1e-8 then true else false;
    end for;
       bool_eps[n] = if (eps_liq[n]<1e-5 and Vdot_l_in<1e-8) then true else false;
    end if;

  /*** MOLE BALANCES ***/
  if n==1 then

     // component balance for vapour
    A*H/n*eps* der(eps_vap[1]*c_v[1,:]) =  upStreamIn.Ndot*actualStream(upStreamIn.x_outflow) + upStreamOut.Ndot*actualStream(upStreamOut.x_outflow) + Ndot_v_transfer[1,:] + Vdot_v_feed[1]*c_v_feed[1,:];
      // component balance for liquid
   A*H/n*eps* der(eps_liq[1]*c_l[1,:]) = downStreamIn.Ndot * actualStream(downStreamIn.x_outflow) + downStreamOut.Ndot* actualStream(downStreamOut.x_outflow) + Ndot_l_transfer[1,:] + Ndot_reac[1,:] + Vdot_l_feed[1]*c_l_feed[1,:] - Vdot_le[1]*c_l[1,:];

    // total mole balance for liquid and vapour
    if eps_liq[1]<1e-5 and Vdot_l_in<1e-8 then
      // bool_eps[1]=true;
       der(eps_liq[1])=0;
   else
     //bool_eps[1]=false;
      A*H/n*eps* der(eps_liq[1]/mediumLiquid[1].v) =  Vdot_l_in*rho_l_in/MM_l_in -  Vdot_l[1]*rho_l[1]/MM_l[1] + sum(Ndot_l_transfer[1,:]) + sum(Ndot_reac[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];
     //    A*H/n*eps* der(eps_liq[1]*rho_l[1]/MM_l[1]) =  Vdot_l_in*rho_l_in/MM_l_in -  Vdot_l[1]*rho_l[1]/MM_l[1] + sum(Ndot_l_transfer[1,:]) + sum(reaction.Ndot[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];
     end if;
         A*H/n*eps* der(eps_vap[1]*rho_v[1]/MM_v[1]) = Vdot_v_in*rho_v_in/MM_v_in -  Vdot_v[1]*rho_v[1]/MM_v[1] + sum(Ndot_v_transfer[1,:]) + Vdot_v_feed[1]*rho_v_feed[1]/MM_v_feed[1] + Ndot_source_startUp[1] + Ndot_source_shutDown[1];
  else

    /** Begin lowest stage (n=1) **/

     // component balance for vapour
    A*H/n*eps* der(eps_vap[1]*c_v[1,:]) = upStreamIn.Ndot*actualStream(upStreamIn.x_outflow) - Vdot_v[1] * c_v[1,:] + Ndot_v_transfer[1,:] + Vdot_v_feed[1]*c_v_feed[1,:] + Ndot_source_startUp[1]*x_v[1,:] + Ndot_source_shutDown[1] * x_v_dummy[:];
      // component balance for liquid
      A*H/n*eps* der(eps_liq[1]*c_l[1,:]) = Vdot_l[2]*c_l[2,:]  +  downStreamOut.Ndot * actualStream(downStreamOut.x_outflow) + Ndot_l_transfer[1,:] + Ndot_reac[1,:] + Vdot_l_feed[1]*c_l_feed[1,:] - Vdot_le[1]*c_l[1,:];

    // total mole balance for liquid and vapour
     if eps_liq[1]<1e-5 and Vdot_l[2]<1e-8 then
       //bool_eps[1]=true;
       der(eps_liq[1])=0;
     else
    // bool_eps[1]=false;
     A*H/n*eps* der(eps_liq[1]/mediumLiquid[1].v) =  Vdot_l[2]*rho_l[2]/MM_l[2] +downStreamOut.Ndot + sum(Ndot_l_transfer[1,:]) + sum(Ndot_reac[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];
     //    A*H/n*eps* der(eps_liq[1]*rho_l[1]/MM_l[1]) =  Vdot_l[2]*rho_l[2]/MM_l[2] -  Vdot_l[1]*rho_l[1]/MM_l[1] + sum(Ndot_l_transfer[1,:]) + sum(reaction.Ndot[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];
       end if;
     A*H/n*eps* der(eps_vap[1]*rho_v[1]/MM_v[1]) = upStreamIn.Ndot -  Vdot_v[1]*rho_v[1]/MM_v[1] + sum(Ndot_v_transfer[1,:]) + Vdot_v_feed[1]*rho_v_feed[1]/MM_v_feed[1] + Ndot_source_startUp[1] + Ndot_source_shutDown[1];
     /** End lowest stage (n=1) **/
    /** Begin stages 2 to n-1 **/
    for j in 2:n-1 loop
      for i in 1:nSV loop
       // component balance for vapour
        A*H/n*eps* der(eps_vap[j]*c_v[j,i]) = Vdot_v[j-1]*c_v[j-1,i] - Vdot_v[j] * c_v[j,i] + Ndot_v_transfer[j,i] + Vdot_v_feed[j]*c_v_feed[j,i] + Ndot_source_startUp[j]*x_v[j,i] + Ndot_source_shutDown[j]*x_v_dummy[i];
      end for;
      for i in 1:nSL loop
        // component balance for liquid
        A*H/n*eps* der(eps_liq[j]*c_l[j,i]) = Vdot_l[j+1]*c_l[j+1,i] - Vdot_l[j] * c_l[j,i] + Ndot_l_transfer[j,i]+ Ndot_reac[j,i]+ Vdot_l_feed[j]*c_l_feed[j,i] + Vdot_le[j-1]*c_l[j-1,i] - Vdot_le[j]*c_l[j,i];
      end for;
     // total mole balance for liquid and vapour
     if eps_liq[j]<1e-5 and Vdot_l[j+1]<1e-8 then
      // bool_eps[j]=true;
      der(eps_liq[j])=0;
     else
       //bool_eps[j]=false;
     A*H/n*eps* der(eps_liq[j]/mediumLiquid[j].v) = Vdot_l[j+1]*rho_l[j+1]/MM_l[j+1] -  Vdot_l[j]*rho_l[j]/MM_l[j]  + sum(Ndot_l_transfer[j,:])+ sum(Ndot_reac[j,:]) + Vdot_l_feed[j]*rho_l_feed[j]/MM_l_feed[j] + Vdot_le[j-1]*rho_l[j-1]/MM_l[j-1] - Vdot_le[j]*rho_l[j]/MM_l[j];
     //   A*H/n*eps* der(eps_liq[j]*rho_l[j]/MM_l[j]) = Vdot_l[j+1]*rho_l[j+1]/MM_l[j+1] -  Vdot_l[j]*rho_l[j]/MM_l[j]  + sum(Ndot_l_transfer[j,:])+ sum(reaction.Ndot[j,:]) + Vdot_l_feed[j]*rho_l_feed[j]/MM_l_feed[j] + Vdot_le[j-1]*rho_l[j-1]/MM_l[j-1] - Vdot_le[j]*rho_l[j]/MM_l[j];
     end if;

   A*H/n*eps* der(eps_vap[j]*rho_v[j]/MM_v[j]) = Vdot_v[j-1]*rho_v[j-1]/MM_v[j-1] -  Vdot_v[j]*rho_v[j]/MM_v[j]  + sum(Ndot_v_transfer[j,:])+ Vdot_v_feed[j]*rho_v_feed[j]/MM_v_feed[j] + Ndot_source_startUp[j] +Ndot_source_shutDown[j];
     end for;
    /** End stages 2 to n-1 **/

    /** Begin highest stage (n=n) **/

     // component balance for vapour
      A*H/n*eps* der(eps_vap[n]*c_v[n,:]) = Vdot_v[n-1]*c_v[n-1,:] + upStreamOut.Ndot*actualStream(upStreamOut.x_outflow) + Ndot_v_transfer[n,:] + Vdot_v_feed[n]*c_v_feed[n,:]  + Ndot_source_startUp[n]*x_v[n,:] + Ndot_source_shutDown[n]*x_v_dummy[:];
      // component balance for liquid
     A*H/n*eps* der(eps_liq[n]*c_l[n,:]) = downStreamIn.Ndot * actualStream(downStreamIn.x_outflow) - Vdot_l[n] * c_l[n,:] + Ndot_l_transfer[n,:]+ Ndot_reac[n,:] + Vdot_l_feed[n]*c_l_feed[n,:] + Vdot_le[n-1]*c_l[n-1,:];
    // total mole balance for liquid and vapour

     if eps_liq[n]<1e-5 and Vdot_l_in<1e-8 then
       //bool_eps[n]=true;
      der(eps_liq[n])=0;
     else
       //bool_eps[n]=false;
     A*H/n*eps* der(eps_liq[n]/mediumLiquid[n].v) = downStreamIn.Ndot -  Vdot_l[n]*rho_l[n]/MM_l[n] + sum(Ndot_l_transfer[n,:])+ sum(Ndot_reac[n,:]) + Vdot_l_feed[n]*rho_l_feed[n]/MM_l_feed[n] + Vdot_le[n-1]*rho_l[n-1]/MM_l[n-1];
     //   A*H/n*eps* der(eps_liq[n]*rho_l[n]/MM_l[n]) = Vdot_l_in*rho_l_in/MM_l_in -  Vdot_l[n]*rho_l[n]/MM_l[n] + sum(Ndot_l_transfer[n,:])+ sum(reaction.Ndot[n,:]) + Vdot_l_feed[n]*rho_l_feed[n]/MM_l_feed[n] + Vdot_le[n-1]*rho_l[n-1]/MM_l[n-1];
     end if;

    A*H/n*eps* der(eps_vap[n]*rho_v[n]/MM_v[n]) = Vdot_v[n-1]*rho_v[n-1]/MM_v[n-1] +upStreamOut.Ndot + sum(Ndot_v_transfer[n,:])+ Vdot_v_feed[n]*rho_v_feed[n]/MM_v_feed[n] + Ndot_source_startUp[n] + Ndot_source_shutDown[n];
  /** End highest stage (n=n) **/
  end if;

    for j in 1:n loop
  //     sum(x_l[j,:])=1;
  //     sum(x_v[j,:])=1;
      eps_vap[j] = 1- eps_liq[j];
    end for;
  /*** ENERGY BALANCE ***/
  if n==1 then
      A*H/n*eps * der(eps_liq[1]*sum(c_l[1,:])*u_l[1]) + A*H/n*(1-eps) * rho_solid[1]*c_solid*der(T[1]) =sum(-Ndot_v_transfer[1,:].* delta_hv[1,:])+ downStreamIn.Ndot*actualStream(downStreamIn.h_outflow) +  downStreamOut.Ndot*actualStream(downStreamOut.h_outflow)  - Qdot_wall[1] + Edot_l_transfer[1]  + Vdot_l_feed[1]*sum(c_l_feed[1,:])*h_l_feed[1] - Vdot_le[1]*sum(c_l[1,:])*h_l[1] +Qdot_reac[1];
      A*H/n*eps * der(eps_vap[1]*sum(c_v[1,:])*u_v[1])  = upStreamIn.Ndot*actualStream(upStreamIn.h_outflow) +  upStreamOut.Ndot*actualStream(upStreamOut.h_outflow) +Edot_v_transfer[1]  + Vdot_v_feed[1]*sum(c_v_feed[1,:])*h_v_feed[1] + Ndot_source_startUp[1]*h_v[1];
  else
    A*H/n*eps * der(eps_liq[1]*sum(c_l[1,:])*u_l[1]) + A*H/n*(1-eps) * rho_solid[1]*c_solid*der(T[1]) = sum(-Ndot_v_transfer[1,:].* delta_hv[1,:])+Vdot_l[2]*sum(c_l[2,:])*h_l[2] +downStreamOut.Ndot*actualStream(downStreamOut.h_outflow) - Qdot_wall[1] + Edot_l_transfer[1]  + Vdot_l_feed[1]*sum(c_l_feed[1,:])*h_l_feed[1] - Vdot_le[1]*sum(c_l[1,:])*h_l[1] + Qdot_reac[1];
      for j in 2:n-1 loop
    A*H/n*eps * der(eps_liq[j]*sum(c_l[j,:])*u_l[j])  + A*H/n*(1-eps) * rho_solid[j]*c_solid*der(T[j]) = sum(-Ndot_v_transfer[j,:].* delta_hv[j,:])+  Vdot_l[j+1]*sum(c_l[j+1,:])*h_l[j+1] - Vdot_l[j]*sum(c_l[j,:])*h_l[j] -Qdot_wall[j]+Edot_l_transfer[j] + Vdot_l_feed[j]*sum(c_l_feed[j,:])*h_l_feed[j] +Vdot_le[j-1]*sum(c_l[j-1,:])*h_l[j-1] - Vdot_le[j]*sum(c_l[j,:])*h_l[j] + Qdot_reac[j];
     end for;
    A*H/n*eps * der(eps_liq[n]*sum(c_l[n,:])*u_l[n])  + A*H/n*(1-eps) * rho_solid[n]*c_solid*der(T[n]) =sum(-Ndot_v_transfer[n,:].* delta_hv[n,:])+ downStreamIn.Ndot*actualStream(downStreamIn.h_outflow) - Vdot_l[n]*sum(c_l[n,:])*h_l[n] - Qdot_wall[n] + Edot_l_transfer[n] + Vdot_l_feed[n]*sum(c_l_feed[n,:])*h_l_feed[n] + Vdot_le[n-1]*sum(c_l[n-1,:])*h_l[n-1]+ Qdot_reac[n];

    A*H/n*eps * der(eps_vap[1]*sum(c_v[1,:])*u_v[1])  =  upStreamIn.Ndot*actualStream(upStreamIn.h_outflow) - Vdot_v[1]*sum(c_v[1,:])*h_v[1] +Edot_v_transfer[1]  + Vdot_v_feed[1]*sum(c_v_feed[1,:])*h_v_feed[1] + Ndot_source_startUp[1]*h_v[1];
    for j in 2:n-1 loop
    A*H/n*eps * der(eps_vap[j]*sum(c_v[j,:])*u_v[j])  =Vdot_v[j-1] * sum(c_v[j-1,:])*mediumVapour[j-1].h  - Vdot_v[j]*sum(c_v[j,:])*h_v[j] +Edot_v_transfer[j] + Vdot_v_feed[j]*sum(c_v_feed[j,:])*h_v_feed[j]  + Ndot_source_startUp[j]*h_v[j];
    end for;
  A*H/n*eps * der(eps_vap[n]*sum(c_v[n,:])*u_v[n]) = Vdot_v[n-1] * sum(c_v[n-1,:])*mediumVapour[n-1].h + upStreamOut.Ndot*actualStream(upStreamOut.h_outflow) +Edot_v_transfer[n] + Vdot_v_feed[n]*sum(c_v_feed[n,:])*h_v_feed[n] + Ndot_source_startUp[n]*h_v[n];
  end if;
    T=T_l;

    //Saturation Pressure
    for j in 1:n loop
      p_sat_bulk[j,:]={mediumLiquid[j].p_sat[i] for i in 1:nSL};
     /// p_sat[j,:]={mediumLiquid[j].p_sat[i] for i in 1:nSL};
        p_sat[j,:]={mediumLiquidStar[j].p_sat[i] for i in 1:nSL};
    end for;

  for j in 1:n loop
      n_mol_L[j] = A*H/n*eps*eps_liq[j]* rho_l_star[j]/MM_l_star[j];
       n_mol_V[j] = A*H/n*eps*eps_vap[j]* mediumVapourStar[j].d/mediumVapourStar[j].MM;
  end for;

    for i in 1:nS-1 loop
      (n_mol_L[1]*x_l_star[1,mapping[i,2]]+n_mol_V[1]*x_v_star[1,mapping[i,1]])/(n_mol_L[1]+n_mol_V[1]) = check[mapping[i,1]];
    end for;

  initial equation




















  equation
    connect(StartUp_signal, startUpReplicator.u) annotation (Line(points={{-97,-61},
            {-87.5,-61},{-73,-61}},      color={255,0,255}));
    connect(endPhase1.u, booleanExpression.y) annotation (Line(points={{-38.8,11},
            {-43.4,11},{-43.4,12},{-49.3,12}}, color={255,0,255}));
    connect(ShutDown_signal, shutDownReplicator.u) annotation (Line(points={{-97,-45},
            {-61.5,-45},{-24.8,-45}}, color={255,0,255}));
      if switchingCondition_Boiling then
        for i in 1:n loop
        internalStepSequence[i].endPhase2=p_bub[i]+friggelfaktor>p_v[i];
        end for;
      else
        for i in 1:n loop
        internalStepSequence[i].endPhase2=true;
        end for;
      end if;
    connect(shutDownReplicator.y, internalStepSequence.shutDown_signal)
      annotation (Line(points={{-15.6,-45},{-6.2,-45},{-6.2,0}}, color={255,0,255}));
    connect(startUpReplicator.y, internalStepSequence.StartUp_signal1)
      annotation (Line(points={{-61.5,-61},{-1.9,-61},{-1.9,0.1}}, color={255,0,255}));
    connect(endPhase1.y, internalStepSequence.endPhase1)
      annotation (Line(points={{-29.6,11},{-24.8,11},{-24.8,11.2},{-19.8,11.2}}, color={255,0,255}));
    annotation (         Icon(graphics),
      Documentation(revisions="<html>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"1\"><tr>
<td><p align=\"center\"><h4>created by </h4></p></td>
<td><p><a href=\"mailto:karin.dietl@tu-harburg.de\">Karin Dietl</a> &AMP; <a href=\"mailto:andreas.joos@tu-harburg.de\">Andreas Joos</a></p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>creation date </h4></p></td>
<td><p>01.01.2009 </p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>revised by </h4></p></td>
<td><p>nobody so far</p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>last revision </h4></p></td>
<td><p>this is an alpha version... </p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>based on </h4></p></td>
<td></td>
</tr>
</table>
<p><br/><br/><br/>Documentation last revised: 18.7.2011</p>
</html>",   info="<html>
<p><h4>General </h4></p>
<p>Stages are counted from the bottom (n=1: lowest stage). The minimum number of stages is n=1. </p>
<p><u><b><font style=\"color: #008000; \">Startup operation</font></b></u></p>
<p>If a rectification column at time t = 0s shall be empty and cold and the starup operation from such an empty and cold state shall be modeled, the boolean parameter &QUOT;considerStartUp&QUOT; has to be set to true (default value = false). An initial pressure has to be provided.</p>
<p>During start-up the inert gas in the column is not modelled. A variable &QUOT;startUp&QUOT; is used in order to determine wether the switching condition on a stage is already fulfilled or not. The switching condition is fulfilled, when the bubble pressure of the mixture attains the initial pressure specified by the user. At this time instant the variable &QUOT;startUp&QUOT; is set to false, vapour is leaving the stage and the equilibrium condition at the phase boundary is valid.</p>
<p>The equations for the liquid phase for start up have to be provided in the extending classes. </p>
<p><u><b><font style=\"color: #008000; \">Medium Models </font></b></u></p>
<p>The liquid medium models and the vapour medium models can differ both in the number of mediums they contain as well as in the substance types. The parameter nSL is the number of substances in the liquid and nSV is the number of substances in the vapour. The parameter nS is the number of substances which are in the liquid as well as in the vapour phase. This parameter has to be supplied by the user. The arrangement of the different substances in the medium models in in theory arbitrary. The parameter mapping has to be used to map the different vectors one to another. </p>
<p>Example: Vapour = {N2, H2O, CO2}, Liquid = {N2, H+, HCO3- H2O, CO2} , mapping = {{1,1},{2,4},{3,5}}.</p>
<p><u><b><font style=\"color: #008000; \">Mole Balances </font></b></u></p>
<p>The mole balances are written separately for vapour and liquid. There exist one mole balance for each component of each stage. The vapour balance is of the following structure: </p>
<p>Mole storage = convective molar flow rate in - convective molar flow rate out + molar flow rate over phase boundary + feed molar flow rate</p>
<p>The liquid balance is of the following structure: </p>
<p>Mole storage = convective molar flow rate in - convective molar flow rate out + molar flow rate over phase boundary + molar flow rate due to reaction + feed molar flow rate</p>
<p><u><b><font style=\"color: #008000; \">Energy Balances </font></b></u></p>
<p>The energy balances are also written separately for vapour and liquid. There exist one energy balance for each stage. The vapour balance is of the following structure: </p>
<p>Energy storage of the vapour = convective enthaply flow rate in - convective enthalpy flow rate out + heat transfer between the phases + enthalpy flow rate from the liquid to the vapour phase - enthalpy flow rate from the vapour to the liquid phase + enthalpy flow rate of the feed</p>
<p>The liquid balance is of the following structure: </p>
<p>Energy storage of the vapour + energy storage of the solid material = convective enthaply flow rate in - convective enthalpy flow rate out + heat transfer to the wall + heat transfer between the phases + enthalpy flow rate from the liquid to the vapour phase - enthalpy flow rate from the vapour to the liquid phase + enthalpy flow rate of the feed</p>
<p><u><b><font style=\"color: #008000; \">Mass Transfer and thermodynamic equilibrium</font></b></u></p>
<p>The mass transfer equations and the equations for the thermodynamic equilibrium are provided in the film model, which is instantiated in the column specific classes <a href=\"Modelica://ThermalSeparation.Components.Columns.StructuredPackedColumn\">StructuredPackedColumn</a>, <a href=\"Modelica://ThermalSeparation.Components.Columns.RandomPackedColumn\">RandomPackedColumn</a>, <a href=\"Modelica://ThermalSeparation.Components.Columns.TrayColumn\">TrayColumn</a> and <a href=\"Modelica://ThermalSeparation.Components.Columns.SprayColumn\">SprayColumn</a>.</p>
</html>"));
  end BaseColumn_newStartUpShutDown;

  model InternalStepSequence

    Modelica_StateGraph2.Step ColdState(
      initialStep=true,
      nOut=1,
      nIn=1,
      use_activePort=false)
      annotation (Placement(transformation(extent={{-36,32},{-28,40}})));
    Modelica_StateGraph2.Transition T1(use_conditionPort=true)
      annotation (Placement(transformation(extent={{-36,8},{-28,16}})));
    Modelica_StateGraph2.Step HeatUpState(
      nOut=1,
      initialStep=false,
      nIn=1) annotation (Placement(transformation(extent={{-4,-4},{4,4}},
          rotation=0,
          origin={-32,-14})));
    Modelica_StateGraph2.Transition T2(use_conditionPort=true)
      annotation (Placement(transformation(extent={{-4,-4},{4,4}},
          rotation=0,
          origin={-32,-38})));
    Modelica_StateGraph2.Step NormalOperation(
      nOut=1,
      initialStep=false,
      nIn=1) annotation (Placement(transformation(
          extent={{-4,-4},{4,4}},
          rotation=90,
          origin={-4,-54})));
    Modelica_StateGraph2.Transition T3(use_conditionPort=true,
      delayedTransition=true,
      waitTime=1)                                              annotation (Placement(
          transformation(
          extent={{-4,-4},{4,4}},
          rotation=180,
          origin={20,-40})));
    Modelica_StateGraph2.Step ShutDown(
      initialStep=false,
      nIn=1,
      nOut=1) annotation (Placement(transformation(
          extent={{-4,-4},{4,4}},
          rotation=180,
          origin={20,12})));
    Modelica_StateGraph2.Transition T4(use_conditionPort=true) annotation (Placement(
          transformation(
          extent={{-4,-4},{4,4}},
          rotation=270,
          origin={-2,60})));
    Modelica.Blocks.Interfaces.BooleanInput shutDown_signal annotation (Placement(
          transformation(
          extent={{-16,-16},{16,16}},
          rotation=90,
          origin={38,-100})));
    Modelica.Blocks.Interfaces.BooleanInput StartUp_signal1 annotation (Placement(
          transformation(
          extent={{-17,-17},{17,17}},
          rotation=90,
          origin={81,-99})));
    Modelica.Blocks.Interfaces.BooleanInput endPhase1 annotation (Placement(
        transformation(
        extent={{-16,-16},{16,16}},
        rotation=0,
        origin={-98,12})));
    Modelica.Blocks.Interfaces.BooleanInput endPhase2 annotation (Placement(
        transformation(
        extent={{-2,-2},{2,2}},
        rotation=0,
        origin={-44,-38})));
  equation
    connect(ColdState.outPort[1], T1.inPort)
      annotation (Line(points={{-32,31.4},{-32,26},{-32,20},{-32,16}},
                                                     color={0,0,0}));
    connect(HeatUpState.inPort[1], T1.outPort)
      annotation (Line(points={{-32,-10},{-32,-10},{-32,7}},color={0,0,0}));
    connect(T2.inPort, HeatUpState.outPort[1])
      annotation (Line(points={{-32,-34},{-32,-32},{-32,-25.4},{-32,-18.6}},
                                                       color={0,0,0}));
    connect(NormalOperation.inPort[1], T2.outPort)
      annotation (Line(points={{-8,-54},{-32,-54},{-32,-43}}, color={0,0,0}));
    connect(NormalOperation.outPort[1], T3.inPort) annotation (Line(points={{0.6,-54},
            {0.6,-54},{20,-54},{20,-44}}, color={0,0,0}));
    connect(ShutDown.inPort[1], T3.outPort)
      annotation (Line(points={{20,8},{20,-22},{20,-35}}, color={0,0,0}));
    connect(T4.outPort, ColdState.inPort[1])
      annotation (Line(points={{-7,60},{-32,60},{-32,40}}, color={0,0,0}));
    connect(T4.inPort, ShutDown.outPort[1])
      annotation (Line(points={{2,60},{20,60},{20,16.6}}, color={0,0,0}));
    connect(shutDown_signal, T3.conditionPort)
      annotation (Line(points={{38,-100},{38,-40},{25,-40}}, color={255,0,255}));
    connect(StartUp_signal1, T4.conditionPort)
      annotation (Line(points={{81,-99},{81,80},{-2,80},{-2,65}}, color={255,0,255}));
  connect(endPhase1, T1.conditionPort)
    annotation (Line(points={{-98,12},{-37,12}}, color={255,0,255}));
  connect(endPhase2, T2.conditionPort)
    annotation (Line(points={{-44,-38},{-44,-38},{-37,-38}},
                                                   color={255,0,255}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-94,94},{94,-94}},
            lineColor={85,170,255},
            lineThickness=1),
          Rectangle(
            extent={{-12,60},{12,40}},
            lineColor={135,135,135},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            radius=60,
            lineThickness=0.5),
          Rectangle(
            extent={{-12,0},{12,-20}},
            lineColor={135,135,135},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            radius=60,
            lineThickness=0.5),
          Rectangle(
            extent={{-12,-60},{12,-80}},
            lineColor={135,135,135},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            radius=60,
            lineThickness=0.5),
          Rectangle(
            extent={{-30,19},{32,24}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid,
            radius=10),
          Rectangle(
            extent={{-30,-41},{32,-36}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid,
            radius=10),
          Line(
            points={{0,40},{0,24}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{0,20},{0,0}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{0,-20},{0,-36},{0,-38}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{0,-40},{0,-56},{0,-60}},
            color={0,0,0},
            thickness=0.5),
          Ellipse(
            extent={{-6,-58},{6,-62}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-6,2},{6,-2}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-6,62},{6,58}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Line(
            points={{0,78},{0,62}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{-94,78},{0,78}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{0,-80},{0,-94}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{-94,50},{-12,50}},
            color={255,0,255},
            thickness=0.5),
          Line(
            points={{-94,-10},{-12,-10}},
            color={255,0,255},
            thickness=0.5),
          Line(
            points={{-94,-70},{-12,-70}},
            color={255,0,255},
            thickness=0.5),
          Line(
            points={{-16,54},{-12,50},{-16,46}},
            color={255,0,255},
            thickness=0.5),
          Line(
            points={{-16,-6},{-12,-10},{-16,-14}},
            color={255,0,255},
            thickness=0.5),
          Line(
            points={{-16,-66},{-12,-70},{-16,-74}},
            color={255,0,255},
            thickness=0.5),
          Text(
            extent={{-84,92},{80,82}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid,
            textString="Step Sequence")}),                         Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end InternalStepSequence;

  partial model FeedColumn_newStartUpShutDown
    "column with optional liquid and or vapour feed inlets"
   //extends BaseColumn;
   extends ThermalSeparation.Components.Columns.BaseClasses.BaseColumn_newStartUpShutDown;
  /***LIQUID FEED ***/

  //protected
  model InternalFeedPort "internal model for conditional stream connectors"
      Interfaces.LiquidPortOut feedLiquidInternal;
      Interfaces.GasPortOut feedVapourInternal;
  end InternalFeedPort;
    InternalFeedPort internalFeedPort[n](feedLiquidInternal(redeclare each
        package   Medium =
            MediumLiquid),feedVapourInternal(redeclare each package Medium =
            MediumVapour));
protected
    parameter Integer[n-numberLiquidFeedsInternal] aux1(each fixed=false);
    parameter Integer[n] aux2(each fixed=false);
    parameter Integer counterL(fixed = false) "counting the stages; equals n at end of algorithm";

public
    parameter Boolean hasLiquidFeed = false "true, if there exist a liquid feed" annotation(Dialog(tab="Feed", group="Liquid Feed"));
    parameter Integer numberLiquidFeeds(min=0,max=n) = 1  annotation(Dialog(enable = hasLiquidFeed,tab="Feed", group="Liquid Feed"));
    parameter Integer[numberLiquidFeeds] stageLiquidFeed={2} "number of stage where feed enters the column"
                                                     annotation(Dialog(enable = hasLiquidFeed, tab="Feed", group="Liquid Feed"));
    final parameter Integer numberLiquidFeedsInternal(min=0,max=n) = if hasLiquidFeed then numberLiquidFeeds else 0;
    Interfaces.LiquidPortIn[numberLiquidFeeds] feedLiquid(redeclare each
      package
        Medium = MediumLiquid) if                        hasLiquidFeed
      annotation (Placement(transformation(extent={{-94,-14},{-74,6}}, rotation=0),
          iconTransformation(extent={{-94,-14},{-74,6}})));

    ThermalSeparation.Utilities.MediumLink mediumLink[n];
    MediumLiquid.BaseProperties mediumLiquidFeed[numberLiquidFeeds](each T0=T_ref,  each p=p_v[n+1],h=actualStream(feedLiquid.h_outflow), x=actualStream(feedLiquid.x_outflow)) if                   hasLiquidFeed;

  /***VAPOUR FEED ***/
  //   ThermalSeparation.Interfaces.GasPortOut[n] feedVapourInternal(redeclare each
  //       package Medium =
  //         MediumVapour,                                                                                       x(start=x_v_start), c(start=c_v));
protected
    parameter Integer[n-numberVapourFeedsInternal] testV(each fixed=false);
    parameter Integer[n] test2V(each fixed=false);
    parameter Integer counterLV(fixed = false) "counting the stages; equals n at end of algorithm";

public
    parameter Boolean hasVapourFeed = false "true, if there exist a liquid feed" annotation(Dialog(tab="Feed", group="Vapour Feed"));
    parameter Integer numberVapourFeeds(min=0,max=n) = 1  annotation(Dialog(enable = hasVapourFeed, tab="Feed", group="Vapour Feed"));
    parameter Integer[numberVapourFeeds] stageVapourFeed={1} "number of stage where feed enters the column"
                                                     annotation(Dialog(enable = hasVapourFeed,tab="Feed", group="Vapour Feed"));
    final parameter Integer numberVapourFeedsInternal(min=0,max=n) = if hasVapourFeed then numberVapourFeeds else 0;
  //   ThermalSeparation.Utilities.LinkVapourSink[n] linkVapour(redeclare each
  //       package Medium =
  //         MediumVapour,   p=p_v[1:n]);
    Interfaces.GasPortIn[numberVapourFeeds] feedVapour(redeclare each package
        Medium = MediumVapour) if                        hasVapourFeed
      annotation (Placement(transformation(extent={{-94,8},{-74,28}}, rotation=0),
          iconTransformation(extent={{-94,8},{-74,28}})));

    ThermalSeparation.Utilities.MediumLink mediumVapourLink[n];
    MediumVapour.BaseProperties mediumVapourFeed[numberVapourFeeds](each T0=T_ref, each p=p_v[n+1],h=actualStream(feedVapour.h_outflow),c=c_v_feed_used, x=actualStream(feedVapour.x_outflow),  x_star=actualStream(feedVapour.x_outflow)) if                   hasVapourFeed;
  //c=c_v_feed,
    ThermalSeparation.Interfaces.GasPortIn[numberVapourFeeds] feedVapour_dummy(redeclare
      each package                                                                                    Medium =
          MediumVapour);
    ThermalSeparation.Interfaces.LiquidPortIn[numberLiquidFeeds] feedLiquid_dummy(redeclare
      each package                                                                                       Medium =
          MediumLiquid);
    ThermalSeparation.Interfaces.GasPortOut[numberVapourFeeds] feedVapour_dummy2(redeclare
      each package                                                                                      Medium =
          MediumVapour) if                                                                                                        not hasVapourFeed;
    ThermalSeparation.Interfaces.LiquidPortOut[numberLiquidFeeds] feedLiquid_dummy2(redeclare
      each package                                                                                         Medium =
          MediumLiquid) if                                                                                                            not hasLiquidFeed;
     SourcesSinks.SourceGas sourceGas[numberVapourFeeds](each use_Flow=false, redeclare
      each package                                                                                  Medium =
          MediumVapour,                                                                                                 each T=293.15,each x = x_v_start_const,each Flow=1) if not hasVapourFeed;
    SourcesSinks.SourceLiquid sourceLiquid[numberLiquidFeeds](each use_Flow=false,redeclare
      each package                                                                                       MediumLiquid =
          MediumLiquid,                                                                                                            each T=293.15,each x=x_l_start_const,each Flow=1) if not hasLiquidFeed;
    SourcesSinks.SinkGas sinkGas[numberVapourFeeds](redeclare each package Medium =
          MediumVapour,                                                                        each p=100000) if not hasVapourFeed;
    SourcesSinks.SinkLiquid sinkLiquid[numberLiquidFeeds](redeclare each
      package                                                                    Medium =
          MediumLiquid,                                                                               each p=100000) if not hasLiquidFeed;

    Real c_v_feed_used[numberVapourFeeds,nSV];


  equation
    /*** link to base class ***/
      for i in 1:numberVapourFeeds loop
      if hasVapourFeed then
        connect(feedVapour_dummy[i], feedVapour[i]);
      else
        connect(feedVapour_dummy[i], sourceGas[i].gasPortOut);
        connect(feedVapour_dummy[i], feedVapour_dummy2[i]);
        connect(feedVapour_dummy2[i],sinkGas[i].gasPortIn);
      end if;
    end for;
    for i in 1:numberLiquidFeeds loop
      if hasLiquidFeed then
        connect(feedLiquid_dummy[i], feedLiquid[i]);
      else
        connect(feedLiquid_dummy[i],sourceLiquid[i].liquidPortOut);
        connect(feedLiquid_dummy[i], feedLiquid_dummy2[i]);
        connect(feedLiquid_dummy2[i],sinkLiquid[i].liquidPortIn);
      end if;
    end for;

  for j in 1:n loop
  Vdot_v_feed[j] = internalFeedPort[j].feedVapourInternal.Ndot/sum(c_v_feed[j,:]);
  c_v_feed[j,:] = inStream(internalFeedPort[j].feedVapourInternal.x_outflow[:])*rho_v_feed[j]/MM_v_feed[j];
  h_v_feed[j] = inStream(internalFeedPort[j].feedVapourInternal.h_outflow);
  rho_v_feed[j] = mediumVapourLink[j].mediumConIn.rho;
  MM_v_feed[j] = mediumVapourLink[j].mediumConIn.MM;
  Vdot_l_feed[j] = internalFeedPort[j].feedLiquidInternal.Ndot/sum(c_l_feed[j,:]);
  c_l_feed[j,:] = inStream(internalFeedPort[j].feedLiquidInternal.x_outflow[:])*rho_l_feed[j]/MM_l_feed[j];
  h_l_feed[j] = inStream(internalFeedPort[j].feedLiquidInternal.h_outflow);
  rho_l_feed[j] = mediumLink[j].mediumConIn.rho;
  MM_l_feed[j] = mediumLink[j].mediumConIn.MM;
  end for;

  /***LIQUID FEED ***/
  //conditional connectors feed
  for j in 1:numberLiquidFeeds loop
    connect(feedLiquid[j],internalFeedPort[stageLiquidFeed[j]].feedLiquidInternal);
    connect(mediumLink[stageLiquidFeed[j]].mediumConIn,mediumLiquidFeed[j].mediumConOut);
  end for;
  for j in 1:n loop
    internalFeedPort[j].feedLiquidInternal.p = p_v[j];
  end for;
  if hasLiquidFeed then
    for j in 1:(n-numberLiquidFeeds) loop
      internalFeedPort[aux1[j]].feedLiquidInternal.h_outflow = h_l[aux1[j]];//1e5;
      internalFeedPort[aux1[j]].feedLiquidInternal.x_outflow = x_l[aux1[j],:];//1/nSL * ones(nSL);//

      mediumLink[aux1[j]].mediumConIn.h = 1;
      mediumLink[aux1[j]].mediumConIn.rho = 1;
      mediumLink[aux1[j]].mediumConIn.MM = 1;
    end for;
    for j in 1:numberLiquidFeeds loop
      internalFeedPort[stageLiquidFeed[j]].feedLiquidInternal.h_outflow = h_l[stageLiquidFeed[j]];//1e5;
      internalFeedPort[stageLiquidFeed[j]].feedLiquidInternal.x_outflow = x_l[stageLiquidFeed[j],:];//1/nSL * ones(nSL);
    end for;
  else
    for j in 1:n loop
      internalFeedPort[j].feedLiquidInternal.h_outflow = h_l[j];//1e5;
      internalFeedPort[j].feedLiquidInternal.x_outflow = x_l[j,:];//nSL * ones(nSL);

      mediumLink[j].mediumConIn.h = 1;
      mediumLink[j].mediumConIn.rho = 1;
      mediumLink[j].mediumConIn.MM = 1;
    end for;
  end if;

  /***VAPOUR FEED ***/
  //conditional connectors feed
  for j in 1:numberVapourFeeds loop
    connect(feedVapour[j],internalFeedPort[stageVapourFeed[j]].feedVapourInternal);
    connect(mediumVapourLink[stageVapourFeed[j]].mediumConIn,mediumVapourFeed[j].mediumConOut);
  end for;
  for j in 1:n loop
    internalFeedPort[j].feedVapourInternal.p = p_v[j];
  end for;
  if hasVapourFeed then
  for j in 1:(n-numberVapourFeeds) loop
    internalFeedPort[testV[j]].feedVapourInternal.h_outflow = h_v[testV[j]];//1e5;//
    internalFeedPort[testV[j]].feedVapourInternal.x_outflow = x_v[testV[j],:];//1/nSV * ones(nSV);//

    mediumVapourLink[testV[j]].mediumConIn.h = 1;
    mediumVapourLink[testV[j]].mediumConIn.rho = 1;
    mediumVapourLink[testV[j]].mediumConIn.MM = 1;
  end for;
    for j in 1:numberVapourFeeds loop
      internalFeedPort[stageVapourFeed[j]].feedVapourInternal.h_outflow = h_v[stageVapourFeed[j]];//1e5;
      internalFeedPort[stageVapourFeed[j]].feedVapourInternal.x_outflow = x_v[stageVapourFeed[j],:];//
    end for;
  else
    for j in 1:n loop
      internalFeedPort[j].feedVapourInternal.h_outflow = h_v[j];//1e5;
      internalFeedPort[j].feedVapourInternal.x_outflow = x_v[j,:];//nSL * ones(nSL);

      mediumVapourLink[j].mediumConIn.h = 1;
      mediumVapourLink[j].mediumConIn.rho = 1;
      mediumVapourLink[j].mediumConIn.MM = 1;
    end for;
  end if;

  for i in 1:numberVapourFeeds loop
      c_v_feed_used[i]= c_v_feed[stageVapourFeed[i],:];
  end for;

  /***LIQUID FEED ***/
  initial algorithm
  aux2 :={i for i in 1:n};
  for j in 1:numberLiquidFeedsInternal loop
      for i in j:n loop
          if stageLiquidFeed[j]==aux2[i] then
              aux2[i]:=0*aux2[i];
          end if;
      end for;
  end for;
  counterL:=1;
  for g in 1:n loop
      if aux2[g]>0 then
          aux1[counterL]:=aux2[g];
          counterL:=counterL + 1;
      end if;
  end for;

  /***VAPOUR FEED ***/
  test2V :={i for i in 1:n};
  for j in 1:numberVapourFeedsInternal loop
      for i in j:n loop
          if stageVapourFeed[j]==test2V[i] then
              test2V[i]:=0*test2V[i];
          end if;
      end for;
  end for;
  counterLV:=1;
  for g in 1:n loop
      if test2V[g]>0 then
          testV[counterLV]:=test2V[g];
          counterLV:=counterLV + 1;
      end if;
  end for;
    annotation (Diagram(graphics),
                         Icon(graphics),
      Documentation(info="<html>
<p>This class provides the equations necessary to describe vapour and / or liquid feeds to the column.</p>
</html>"));
  end FeedColumn_newStartUpShutDown;
end BaseClasses;
