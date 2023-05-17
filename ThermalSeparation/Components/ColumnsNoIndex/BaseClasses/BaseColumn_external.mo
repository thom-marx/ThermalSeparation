within ThermalSeparation.Components.ColumnsNoIndex.BaseClasses;
partial model BaseColumn_external "Gesamtmolbilanz, x, h, Ndot im Konnektor"
  import ThermalSeparation;

outer ThermalSeparation.SystemTS systemTS;
parameter SI.Temperature T_ref = systemTS.T_ref "reference temperature" annotation(Dialog(tab="Advanced"));

  parameter Integer n(min=1)=1
    "packed column: number of discrete elements in the section; plate column: number of trays in one section";

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
  final parameter SI.MoleFraction x_l_start[n,nSL](fixed=false);
  final parameter SI.MoleFraction x_v_start[n,nSV](fixed=false);
  parameter Real x_total_start[nSV]=fill(1/nSV,nSV)
    "total mole fraction in system (vapour and liquid), component ordering as in vapour medium"
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

//Result-Record
Results results(n=n, nSL=nSL, nSV=nSV, T_l=T_l, x_l=x_l, x_v=x_v,  x_l_star=x_l_star, x_v_star=x_v_star,c_l = c_l, c_v=c_v, p_v=p_v, Vdot_l=Vdot_l, Vdot_v=Vdot_v, eps_liq = eps_liq, startUp=startUp);

//Equilibrium or Non-Equilibrium model?
  parameter Boolean EQ= false
    "equilibrium model is used, no mass transfer, value provided by film model" annotation(Dialog(enable=false));
//Medium
parameter Integer mapping[nS,2] = {{i,i} for i in 1:nS}
    "parameter to map the different medium vectors one to another";
parameter Boolean inertVapour[nSV] = fill(false,nSV)
    "true for each component which is inert in the vapour phase";
parameter Boolean inertLiquid[nSL] = fill(false,nSL)
    "true for each component which is inert in the liquid phase";
replaceable package MediumVapour =
    ThermalSeparation.Media.BaseMediumVapour
    "medium to be used in vapour phase"                                                                                                  annotation(choicesAllMatching);
 MediumVapour.BaseProperties mediumVapour[n](c=c_v, each T0=T_ref, p=p_v[1:n],  T=T_v,  x=x_v,   x_star=x_v_star);
 MediumVapour.BaseProperties mediumVapourIn(c=c_v_in, T0=T_ref, p=upStreamIn.p, T=T_v_in, x=x_v_in, x_star=x_v_in);
  MediumVapour.BaseProperties mediumVapourStar[n]( c=c_v,  each T0=T_ref,  p=p_v[1:n],  T=T_star, x=x_v_star, x_star=x_v_star);

replaceable package MediumLiquid =
    ThermalSeparation.Media.BaseMediumLiquid
    "medium to be used in liquid phase"                                                                                                  annotation(choicesAllMatching);
       MediumLiquid.BaseProperties mediumLiquid[n](each T0=T_ref, p=p_hyd[1:n], T=T_l, x= x_l, h=h_l);
 MediumLiquid.BaseProperties mediumLiquidIn(T0=T_ref,  p=p_hyd[n+1], T=T_l_in, x=x_l_in,h=h_l_in);
 MediumLiquid.BaseProperties mediumLiquidStar[n](each T0=T_ref, p=p_hyd[1:n], T=T_star, x=x_l_star,h=h_l_star);
 MediumLiquid.ActivityCoefficient activityCoeff[n]( T=T_star, x_l=x_l_star);
 MediumVapour.EvaporationEnthalpy evapEnthalpy[n]( p=p_hyd[1:n],  T=T_v);
 constant Boolean h_evap_medium = MediumVapour.delta_hv_medium;
 ThermalSeparation.Units.MolarEnthalpy delta_hv[n,nSV] = if h_evap_medium then zeros(n,nSV) else evapEnthalpy.h;
  Real h_l_star[n];
  parameter Integer nS(min=2)
    "number of species which are equal in vapour and liquid phase";
  final parameter Integer nL=MediumLiquid.nSubstance -nS
    "number of additional substances which are only in liquid phase";
  final parameter Integer nV = MediumVapour.nSubstance -nS
    "number of additional substances which are only in the vapour phase";
  final parameter Integer nSL = MediumLiquid.nSubstance;
  final parameter Integer nSV = MediumVapour.nSubstance;

  /*** vapour properties ***/
  SI.Density rho_v[n]= mediumVapour.d "mixture vapour density";
  SI.Density rho_v_in =  mediumVapourIn.d;
  SI.MolarMass MM_v[n]( start=0.028*ones(n))= mediumVapour.MM
    "molar mass of the vapour mixture ";
  SI.MolarMass MM_v_in( start=0.03) = mediumVapourIn.MM;
  ThermalSeparation.Units.MolarEnthalpy h_v[n] = mediumVapour.h;
  ThermalSeparation.Units.MolarEnthalpy h_v_in = mediumVapourIn.h;
  SI.MolarInternalEnergy u_v[n](stateSelect=StateSelect.prefer)= mediumVapour.u;
  SI.Concentration c_v_star[n,nSV];
  SI.Density rho_v_star[n] = mediumVapourStar.d;

  /*** liquid properties ***/
  SI.Density rho_l[n] = mediumLiquid.d "mixture liquid density";
  SI.Density rho_l_in = mediumLiquidIn.d;
  SI.Density rho_l_star[n]= mediumLiquidStar.d;
  SI.MolarMass MM_l[n](start=fill(0.018,n))= mediumLiquid.MM
    "molar mass of the liquid mixture";
  SI.MolarMass MM_l_in= mediumLiquidIn.MM;
  SI.MolarMass MM_l_star[n]= mediumLiquidStar.MM;
  ThermalSeparation.Units.MolarEnthalpy h_l[n]= mediumLiquid.h;
  ThermalSeparation.Units.MolarEnthalpy h_l_in=mediumLiquidIn.h;
  SI.MolarInternalEnergy u_l[n](stateSelect=StateSelect.prefer) =  mediumLiquid.u;

//Variables upStream
  SI.Concentration c_v_in[nSV];
  SI.Concentration c_v[n,nSV](stateSelect=StateSelect.default) annotation(Dialog(group="Initialization",showStartAttribute=true));
  SI.MoleFraction x_v_in[nSV];
  SI.MoleFraction x_v[n,nSV](start=x_v_start);
  SI.VolumeFlowRate Vdot_v_in(nominal=1e-2);
  SI.VolumeFlowRate Vdot_v[n](nominal=fill(1e-2,n));
  SI.Temperature T_v_in;

  SI.Pressure p_v[n+1]( start=fill(1.0e5,n+1))
    "p_v[j] = pressure on the j-th stage, p_v[n+1] is the pressure in the first element of the sucesseding component";
  SI.Temperature T_v[n](nominal=fill(350,n),start=fill(350,n));

  //Variables downStream
  SI.Concentration c_l_in[nSL]
    "molar concentration in the liquid at the liquid outlet of each stage";
  SI.Concentration c_l[n,nSL](stateSelect=StateSelect.default) annotation(Dialog(group="Initialization",showStartAttribute=true));

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
  ThermalSeparation.Interfaces.GasPortIn          upStreamIn(
                                                    redeclare package Medium =
        MediumVapour)               annotation (Placement(transformation(extent={{-80,
            -100},{-60,-80}},      rotation=0), iconTransformation(extent={{-80,
            -100},{-60,-80}})));
  ThermalSeparation.Interfaces.GasPortOut          upStreamOut(
                                                      redeclare package Medium =
        MediumVapour)                 annotation (Placement(transformation(
          extent={{-80,80},{-60,100}}, rotation=0), iconTransformation(extent={
            {-80,80},{-60,100}})));
  ThermalSeparation.Interfaces.LiquidPortIn          downStreamIn(
                                       redeclare package Medium=MediumLiquid)
                                       annotation (Placement(transformation(
          extent={{60,80},{80,100}}, rotation=0), iconTransformation(extent={{
            60,80},{80,100}})));
  ThermalSeparation.Interfaces.LiquidPortOut          downStreamOut(
                                         redeclare package Medium =
        MediumLiquid)                    annotation (Placement(transformation(
          extent={{60,-100},{80,-80}}, rotation=0), iconTransformation(extent={
            {60,-100},{80,-80}})));

  //initial equation for eps_liq is supplied in the extending class!
  SI.VolumeFraction eps_liq[        n]( stateSelect=StateSelect.default)
    "liquid volume fraction";
  SI.VolumeFraction eps_vap[        n](start=fill(0.99,n))
    "vapour volume fraction";
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

/*** Variables needed for FluidProp interface ***/
Real n_tot[n,nS];
SI.MoleFraction x_vap_liq[n,nS] "total molar fractions";

/*** thermodynamic equilibrium ***/
 replaceable model ThermoEquilibrium =
      ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid                                  constrainedby ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium
    "model for phase equilibrium"                                                         annotation(choicesAllMatching=true);
    ThermoEquilibrium bubblePressure[n](x_vap_liq=x_vap_liq,each nS=nS,  each mapping =                            mapping,
    redeclare replaceable package MediumVapour =   MediumVapour,
  redeclare replaceable package MediumLiquid =
     MediumLiquid, p=ones(n)*p_initial, T=T_l, x_v=x_v, x_l=x_l, p_sat=p_sat_bulk,  v_v=MM_v./rho_v);

 // replaceable model ThermoEquilibrium =
 //      ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid                                  constrainedby
 //    ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium
 //    "model for phase equilibrium"                                                         annotation(choicesAllMatching=true);
 //    ThermoEquilibrium bubblePressure[n](x_vap_liq=x_vap_liq, each nS=nS,  each mapping =                            mapping,
 //    redeclare replaceable package MediumVapour =   MediumVapour,
 //  redeclare replaceable package MediumLiquid =
 //     MediumLiquid, p=ones(n)*p_initial, T=T_l, x_v=x_v, x_l=x_l, p_sat=p_sat_bulk,  v_v=MM_v./rho_v);

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
  parameter Boolean considerStartUp = false
    "true if StartUp is to be considered" annotation(Dialog(tab="StartUp"));

  SI.Pressure p_initial = 1e5;//upStreamOut.p "initial pressure in column";
  parameter Real friggelfaktor = 0.0002e5 annotation(Dialog(tab="StartUp"));
  parameter Real k=0.2e-3 "large value for steep omega" annotation(Dialog(tab="StartUp"));
  SI.Pressure p_bub[n]= bubblePressure.p_bubble "mixture bubble pressure";
  SI.Pressure p_hyd[n+1] "hydraulic pressure";
  Real omega[n];
  Boolean startUp[n](start=fill(true,n),fixed=false);
  Real Ndot_source_startUp[n]
    "dummy molar flow rate to account for discharge of inert gas during startUp";

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

SI.MolarMass MM_v_star[n](start=fill(0.02,n))=mediumVapourStar.MM;
Real check[nS-1];

//         ThermalSeparation.Utilities.LimPID_Input PID[
//                                               n](each u_s = p_initial,u_m = if considerStartUp then p_v[1:n] else fill(p_initial,n),
//       yMax=100,
//       k=0.005,
//      initType=Modelica.Blocks.Types.Init.InitialOutput)   annotation (Placement(transformation(extent={{-24,-12},{-4,8}})));

        ThermalSeparation.Utilities.LimPID_Input PID[
                                              n](each u_s = p_initial,
      each yMax=100,
      each k=0.5,
     each initType=Modelica.Blocks.Types.Init.InitialOutput)   annotation (Placement(transformation(extent={{-24,-12},{-4,8}})));

/*** changes to reduce index ***/
Modelica.Units.SI.AmountOfSubstance n_mol[n,nS];
Modelica.Units.SI.AmountOfSubstance n_mol_L[n](stateSelect=StateSelect.prefer);
Modelica.Units.SI.AmountOfSubstance n_mol_V[n](stateSelect=StateSelect.prefer);
Real n_mol_L_i[n,nSL](stateSelect=StateSelect.prefer);
Real n_mol_V_i[n,nSV](stateSelect=StateSelect.prefer);

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
            PID[j].u_m = if considerStartUp and startUp[j] then p_v[j] else if considerStartUp and not startUp[j] then p_initial else p_initial;

      end for;
  /***monitoring ***/
  for j in 1:n loop
  X_v[j,:]=x_v[j,:]/MM_v[j].*MediumVapour.MMX[:];
  X_l[j,:]=x_l[j,:]/MM_l[j].*MediumLiquid.MMX[:];
    n_i_liq[j,:]=c_l[j,:].*eps_liq[j]*A*H/n*eps;
    n_i_vap[j,:]=c_v[j,:].*eps_vap[j]*A*H/n*eps;
  end for;

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
  upStreamIn.Ndot = Ndot_v_in;
  upStreamOut.Ndot = -Ndot_v[n];
 //upStreamOut.p_medium = p_hyd[n];
  //downstream
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
    p_hyd[n+1] = p_v[n+1];

else  startUp = fill(false,n);
      omega = ones(n);
      p_hyd = p_v;
      Ndot_source_startUp=zeros(n);
end if;

  for j in 1:n loop
  x_vap_liq[j,:] = n_tot[j,:]./sum(n_tot[j,:]);
    for i in 1:nS loop
    n_tot[j,i]= (c_v[j,mapping[i,1]].* (1-eps_liq[j]) + c_l[j,mapping[i,2]].* eps_liq[j])*H*A/n;
    end for;
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
/*** MOLE BALANCES ***/
if n==1 then

   // component balance for vapour
  der(n_mol_V_i[1,:]) =  upStreamIn.Ndot*actualStream(upStreamIn.x_outflow) + upStreamOut.Ndot*actualStream(upStreamOut.x_outflow) + Ndot_v_transfer[1,:] + Vdot_v_feed[1]*c_v_feed[1,:] + Ndot_source_startUp[1] * x_v[1,:];
    // component balance for liquid
 der(n_mol_L_i[1,:]) = downStreamIn.Ndot * actualStream(downStreamIn.x_outflow) + downStreamOut.Ndot* actualStream(downStreamOut.x_outflow) + Ndot_l_transfer[1,:] + Ndot_reac[1,:] + Vdot_l_feed[1]*c_l_feed[1,:] - Vdot_le[1]*c_l[1,:];

  // total mole balance for liquid and vapour
   bool_eps[1]=false;
   der(n_mol_L[1]) =  Vdot_l_in*rho_l_in/MM_l_in -  Vdot_l[1]*rho_l[1]/MM_l[1] + sum(Ndot_l_transfer[1,:]) + sum(Ndot_reac[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];
   //    A*H/n*eps* der(eps_liq[1]*rho_l[1]/MM_l[1]) =  Vdot_l_in*rho_l_in/MM_l_in -  Vdot_l[1]*rho_l[1]/MM_l[1] + sum(Ndot_l_transfer[1,:]) + sum(reaction.Ndot[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];

      der(n_mol_V[1]) = Vdot_v_in*rho_v_in/MM_v_in -  Vdot_v[1]*rho_v[1]/MM_v[1] + sum(Ndot_v_transfer[1,:]) + Vdot_v_feed[1]*rho_v_feed[1]/MM_v_feed[1] + Ndot_source_startUp[1];
else

  /** Begin lowest stage (n=1) **/

   // component balance for vapour
  der(n_mol_V_i[1,:]) = upStreamIn.Ndot*actualStream(upStreamIn.x_outflow) - Vdot_v[1] * c_v[1,:] + Ndot_v_transfer[1,:] + Vdot_v_feed[1]*c_v_feed[1,:] + Ndot_source_startUp[1]*x_v[1,:];
    // component balance for liquid
   der(n_mol_L_i[1,:]) = Vdot_l[2]*c_l[2,:]  +  downStreamOut.Ndot * actualStream(downStreamOut.x_outflow) + Ndot_l_transfer[1,:] + Ndot_reac[1,:] + Vdot_l_feed[1]*c_l_feed[1,:] - Vdot_le[1]*c_l[1,:];

  // total mole balance for liquid and vapour

   bool_eps[1]=false;
   der(n_mol_L[1]) =  Vdot_l[2]*rho_l[2]/MM_l[2] +downStreamOut.Ndot + sum(Ndot_l_transfer[1,:]) + sum(Ndot_reac[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];
   //    A*H/n*eps* der(eps_liq[1]*rho_l[1]/MM_l[1]) =  Vdot_l[2]*rho_l[2]/MM_l[2] -  Vdot_l[1]*rho_l[1]/MM_l[1] + sum(Ndot_l_transfer[1,:]) + sum(reaction.Ndot[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];

   der(n_mol_V[1]) = upStreamIn.Ndot -  Vdot_v[1]*rho_v[1]/MM_v[1] + sum(Ndot_v_transfer[1,:]) + Vdot_v_feed[1]*rho_v_feed[1]/MM_v_feed[1] + Ndot_source_startUp[1];
  /** End lowest stage (n=1) **/

  /** Begin stages 2 to n-1 **/
  for j in 2:n-1 loop
    for i in 1:nSV loop
     // component balance for vapour
     der(n_mol_V_i[j,i]) = Vdot_v[j-1]*c_v[j-1,i] - Vdot_v[j] * c_v[j,i] + Ndot_v_transfer[j,i] + Vdot_v_feed[j]*c_v_feed[j,i] + Ndot_source_startUp[j]*x_v[j,i];
    end for;
    for i in 1:nSL loop
      // component balance for liquid
     der(n_mol_L_i[j,i]) = Vdot_l[j+1]*c_l[j+1,i] - Vdot_l[j] * c_l[j,i] + Ndot_l_transfer[j,i]+ Ndot_reac[j,i]+ Vdot_l_feed[j]*c_l_feed[j,i] + Vdot_le[j-1]*c_l[j-1,i] - Vdot_le[j]*c_l[j,i];
    end for;
   // total mole balance for liquid and vapour
     bool_eps[j]=false;
   der(n_mol_L[j]) = Vdot_l[j+1]*rho_l[j+1]/MM_l[j+1] -  Vdot_l[j]*rho_l[j]/MM_l[j]  + sum(Ndot_l_transfer[j,:])+ sum(Ndot_reac[j,:]) + Vdot_l_feed[j]*rho_l_feed[j]/MM_l_feed[j] + Vdot_le[j-1]*rho_l[j-1]/MM_l[j-1] - Vdot_le[j]*rho_l[j]/MM_l[j];
   //   A*H/n*eps* der(eps_liq[j]*rho_l[j]/MM_l[j]) = Vdot_l[j+1]*rho_l[j+1]/MM_l[j+1] -  Vdot_l[j]*rho_l[j]/MM_l[j]  + sum(Ndot_l_transfer[j,:])+ sum(reaction.Ndot[j,:]) + Vdot_l_feed[j]*rho_l_feed[j]/MM_l_feed[j] + Vdot_le[j-1]*rho_l[j-1]/MM_l[j-1] - Vdot_le[j]*rho_l[j]/MM_l[j];

  der(n_mol_V[j]) = Vdot_v[j-1]*rho_v[j-1]/MM_v[j-1] -  Vdot_v[j]*rho_v[j]/MM_v[j]  + sum(Ndot_v_transfer[j,:])+ Vdot_v_feed[j]*rho_v_feed[j]/MM_v_feed[j] + Ndot_source_startUp[j];
   end for;
  /** End stages 2 to n-1 **/

  /** Begin highest stage (n=n) **/

   // component balance for vapour

   der(n_mol_V_i[n,:]) = Vdot_v[n-1]*c_v[n-1,:] + upStreamOut.Ndot*actualStream(upStreamOut.x_outflow) + Ndot_v_transfer[n,:] + Vdot_v_feed[n]*c_v_feed[n,:]  + Ndot_source_startUp[n]*x_v[n,:];

    // component balance for liquid
   der(n_mol_L_i[n,:]) = downStreamIn.Ndot * actualStream(downStreamIn.x_outflow) - Vdot_l[n] * c_l[n,:] + Ndot_l_transfer[n,:]+ Ndot_reac[n,:] + Vdot_l_feed[n]*c_l_feed[n,:] + Vdot_le[n-1]*c_l[n-1,:];

  // total mole balance for liquid and vapour

     bool_eps[n]=false;
   der(n_mol_L[n]) = downStreamIn.Ndot -  Vdot_l[n]*rho_l[n]/MM_l[n] + sum(Ndot_l_transfer[n,:])+ sum(Ndot_reac[n,:]) + Vdot_l_feed[n]*rho_l_feed[n]/MM_l_feed[n] + Vdot_le[n-1]*rho_l[n-1]/MM_l[n-1];
   //   A*H/n*eps* der(eps_liq[n]*rho_l[n]/MM_l[n]) = Vdot_l_in*rho_l_in/MM_l_in -  Vdot_l[n]*rho_l[n]/MM_l[n] + sum(Ndot_l_transfer[n,:])+ sum(reaction.Ndot[n,:]) + Vdot_l_feed[n]*rho_l_feed[n]/MM_l_feed[n] + Vdot_le[n-1]*rho_l[n-1]/MM_l[n-1];

   der(n_mol_V[n]) = Vdot_v[n-1]*rho_v[n-1]/MM_v[n-1] +upStreamOut.Ndot + sum(Ndot_v_transfer[n,:])+ Vdot_v_feed[n]*rho_v_feed[n]/MM_v_feed[n] + Ndot_source_startUp[n];
/** End highest stage (n=n) **/
end if;

  for j in 1:n loop
//     sum(x_l[j,:])=1;
//     sum(x_v[j,:])=1;
    eps_vap[j] = 1- eps_liq[j];
  end for;
/*** ENERGY BALANCE ***/
if n==1 then
    der(n_mol_L[1]*u_l[1]) =sum(-Ndot_v_transfer[1,:].* delta_hv[1,:])+ downStreamIn.Ndot*actualStream(downStreamIn.h_outflow) +  downStreamOut.Ndot*actualStream(downStreamOut.h_outflow)  - Qdot_wall[1] + Edot_l_transfer[1]  + Vdot_l_feed[1]*sum(c_l_feed[1,:])*h_l_feed[1] - Vdot_le[1]*sum(c_l[1,:])*h_l[1] +Qdot_reac[1];
    der(n_mol_V[1]*u_v[1])  = upStreamIn.Ndot*actualStream(upStreamIn.h_outflow) +  upStreamOut.Ndot*actualStream(upStreamOut.h_outflow) +Edot_v_transfer[1]  + Vdot_v_feed[1]*sum(c_v_feed[1,:])*h_v_feed[1] + Ndot_source_startUp[1]*h_v[1];
else
  der(n_mol_L[1]*u_l[1]) = sum(-Ndot_v_transfer[1,:].* delta_hv[1,:])+Vdot_l[2]*sum(c_l[2,:])*h_l[2] +downStreamOut.Ndot*actualStream(downStreamOut.h_outflow) - Qdot_wall[1] + Edot_l_transfer[1]  + Vdot_l_feed[1]*sum(c_l_feed[1,:])*h_l_feed[1] - Vdot_le[1]*sum(c_l[1,:])*h_l[1] + Qdot_reac[1];
    for j in 2:n-1 loop
  der(n_mol_L[j]*u_l[j]) = sum(-Ndot_v_transfer[j,:].* delta_hv[j,:])+  Vdot_l[j+1]*sum(c_l[j+1,:])*h_l[j+1] - Vdot_l[j]*sum(c_l[j,:])*h_l[j] -Qdot_wall[j]+Edot_l_transfer[j] + Vdot_l_feed[j]*sum(c_l_feed[j,:])*h_l_feed[j] +Vdot_le[j-1]*sum(c_l[j-1,:])*h_l[j-1] - Vdot_le[j]*sum(c_l[j,:])*h_l[j] + Qdot_reac[j];
   end for;
  der(n_mol_L[n]*u_l[n])=sum(-Ndot_v_transfer[n,:].* delta_hv[n,:])+ downStreamIn.Ndot*actualStream(downStreamIn.h_outflow) - Vdot_l[n]*sum(c_l[n,:])*h_l[n] - Qdot_wall[n] + Edot_l_transfer[n] + Vdot_l_feed[n]*sum(c_l_feed[n,:])*h_l_feed[n] + Vdot_le[n-1]*sum(c_l[n-1,:])*h_l[n-1]+ Qdot_reac[n];

  der(n_mol_V[1]*u_v[1])  =  upStreamIn.Ndot*actualStream(upStreamIn.h_outflow) - Vdot_v[1]*sum(c_v[1,:])*h_v[1] +Edot_v_transfer[1]  + Vdot_v_feed[1]*sum(c_v_feed[1,:])*h_v_feed[1] + Ndot_source_startUp[1]*h_v[1];
  for j in 2:n-1 loop
  der(n_mol_V[j]*u_v[j])  =Vdot_v[j-1] * sum(c_v[j-1,:])*mediumVapour[j-1].h  - Vdot_v[j]*sum(c_v[j,:])*h_v[j] +Edot_v_transfer[j] + Vdot_v_feed[j]*sum(c_v_feed[j,:])*h_v_feed[j]  + Ndot_source_startUp[j]*h_v[j];
  end for;
  der(n_mol_V[n]*u_v[n]) = Vdot_v[n-1] * sum(c_v[n-1,:])*mediumVapour[n-1].h + upStreamOut.Ndot*actualStream(upStreamOut.h_outflow) +Edot_v_transfer[n] + Vdot_v_feed[n]*sum(c_v_feed[n,:])*h_v_feed[n] + Ndot_source_startUp[n]*h_v[n];
end if;
  T=T_l;

  //Saturation Pressure
  for j in 1:n loop
    p_sat_bulk[j,:]={mediumLiquid[j].p_sat[i] for i in 1:nSL};
   /// p_sat[j,:]={mediumLiquid[j].p_sat[i] for i in 1:nSL};
      p_sat[j,:]={mediumLiquidStar[j].p_sat[i] for i in 1:nSL};
  end for;

for j in 1:n loop
    n_mol_V[j] = A*H/n*eps*eps_vap[j]* sum(c_v[j,:]);
    n_mol_L[j] = A*H/n*eps*eps_liq[j]* sum(c_l[j,:]);
    n_mol_L_i[j,:]=n_mol_L[j]*x_l[j,:];
    n_mol_V_i[j,:]=n_mol_V[j]*x_v[j,:];
    for i in 1:nS loop
    n_mol[j,i]=A*H/n*eps*eps_vap[j]* c_v[j,mapping[i,1]]+A*H/n*eps*eps_liq[j]*c_l[j,mapping[i,1]];
    end for;
end for;

  for i in 1:nS-1 loop
    (n_mol_L[1]*x_l_star[1,mapping[i,2]]+n_mol_V[1]*x_v_star[1,mapping[i,1]])/(n_mol_L[1]+n_mol_V[1]) = check[mapping[i,1]];
  end for;

initial equation

  annotation (Diagram(graphics),
                       Icon(graphics),
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
</html>", info="<html>
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
end BaseColumn_external;
