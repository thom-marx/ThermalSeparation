within ThermalSeparation.BalanceEquations.Base.Equilibrium;
partial model BaseBalanceEquationsEq

extends ThermalSeparation.BalanceEquations.Base.BaseBalanceEquations;

  parameter Integer n(min=1)
    "packed column: number of discrete elements in the section; plate column: number of trays in one section";
  parameter Integer nS
    "number of species which are equal in vapour and liquid phase";
  parameter Integer nSL=MediumLiquid.nSubstance;
  parameter Integer nSV=MediumVapour.nSubstance;
  parameter Integer mapping[nS,2] = {{i,i} for i in 1:nS}
    "parameter to map the different medium vectors one to another";

  replaceable package MediumLiquid = Media.BaseMediumLiquid annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  replaceable package MediumVapour = Media.BaseMediumVapour annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

  parameter Boolean inertVapour[nSV] = fill(false,nSV)
    "true for each component which is inert in the vapour phase";
  parameter Boolean inertLiquid[nSL] = fill(false,nSL)
    "true for each component which is inert in the liquid phase";

  /*** interface medium models and properties ***/
  SI.Temperature T_star[n];
  input SI.MoleFraction x_l_star[n,nSL]; // is this really an input or does it come from the film model?
  output SI.MoleFraction c_l_star[n,nSL];
  input SI.MoleFraction x_v_star[n,nSV]; // is this really an input or does it come from the film model?
  output SI.MoleFraction c_v_star[n,nSV];

  /*** vapour properties ***/
  input MediumVapour.ThermodynamicProperties[n] propsVap;
  input MediumVapour.ThermodynamicProperties propsVapIn;
  input MediumVapour.ThermodynamicState[n] stateVap;
   SI.Pressure p_sat[n,nSL];
   SI.Density rho_v[n]=propsVap.d "mixture vapour density";
   SI.Density rho_v_in=propsVapIn.d;
   SI.MolarMass MM_v[n](start=0.028*ones(n))=propsVap.MM
    "molar mass of the vapour mixture ";
   SI.MolarMass MM_v_in=propsVapIn.MM;
   ThermalSeparation.Units.MolarEnthalpy h_v[n]=propsVap.h;
   ThermalSeparation.Units.MolarEnthalpy h_v_in=propsVapIn.h;
  input SI.Pressure p_v[n+1];
  input SI.Pressure p_sat_bulk[n,nSL];

  /*** liquid properties ***/
  input MediumLiquid.ThermodynamicProperties[n] propsLiq;
  input MediumLiquid.ThermodynamicProperties propsLiqIn;
  input MediumLiquid.ThermodynamicState[n] stateLiq;
   SI.Density rho_l[n]=propsLiq.d "mixture liquid density";
   SI.Density rho_l_in=propsLiqIn.d;
   SI.MolarMass MM_l[n](start=fill(0.018,n))=propsLiq.MM
    "molar mass of the liquid mixture";
   SI.MolarMass MM_l_in=propsLiqIn.MM;
   ThermalSeparation.Units.MolarEnthalpy h_l[n]=propsLiq.h;
   ThermalSeparation.Units.MolarEnthalpy h_l_in=propsLiqIn.h;
   parameter Modelica.SIunits.Temperature T_ref;

  /*** variables upStream ***/
  input SI.Concentration c_v[n,nSV](stateSelect=StateSelect.default);
   SI.Concentration c_v_in[nSV]=propsVapIn.c;
   SI.MoleFraction x_v_in[nSV]=propsVapIn.x;
   // SI.MoleFraction x_v[n,nSV]=propsVap.x;
   SI.MoleFraction x_v[n,nSV];
  input SI.VolumeFlowRate Vdot_v_in(nominal=1e-2);
  input SI.VolumeFlowRate Vdot_v[n](nominal=fill(1e-2,n));
  input SI.MoleFraction x_upStreamIn_act[nSV];
  input SI.MoleFraction x_upStreamOut_act[nSV];
  input ThermalSeparation.Units.MolarEnthalpy h_upStreamIn_act;
  input ThermalSeparation.Units.MolarEnthalpy h_upStreamOut_act;
  input SI.MoleFraction x_vap_liq[n,nS];

  /*** variables downStream ***/
  input SI.Concentration c_l[n,nSL](stateSelect=StateSelect.default);
   SI.Concentration c_l_in[nSL]=propsLiqIn.c
    "molar concentration in the liquid at the liquid outlet of each stage";
   SI.MoleFraction x_l_in[nSL]=propsLiqIn.x;
   SI.MoleFraction x_l[n,nSL]=propsLiq.x;
  input SI.VolumeFlowRate Vdot_l_in(nominal=1e-4);
  input SI.VolumeFlowRate Vdot_l[n](nominal=fill(1e-4,n));
  input SI.MoleFraction x_downStreamIn_act[nSL];
  input SI.MoleFraction x_downStreamOut_act[nSL];
  input ThermalSeparation.Units.MolarEnthalpy h_downStreamIn_act;
  input ThermalSeparation.Units.MolarEnthalpy h_downStreamOut_act;

  /*** feed variables ***/
  input SI.VolumeFlowRate Vdot_v_feed[n];
  input SI.Concentration c_v_feed[n,nSV];
  input SI.SpecificEnthalpy h_v_feed[n];
  input SI.VolumeFlowRate Vdot_l_feed[n];
  input SI.Concentration c_l_feed[n,nSL];
  input SI.SpecificEnthalpy h_l_feed[n];
  input SI.Density rho_l_feed[n];
  input SI.Density rho_v_feed[n];
  input SI.MolarMass MM_l_feed[n];
  input SI.MolarMass MM_v_feed[n];

        SI.MolarFlowRate Ndot_v_transfer[n,nSV](start=fill(-0.1,n,nSV)); // calculated in film model/balance equations
  input SI.MolarFlowRate Ndot_l_transfer[n,nSL](start=fill(0.1,n,nSL)); // calculated in film model/balance equations
        SI.HeatFlowRate Edot_l_transfer[n]; // calculated in film model/balance equations
        SI.HeatFlowRate Edot_v_transfer[n]; // calculated in film model/balance equations
  input SI.MolarFlowRate Ndot_reac[n,nSL];
  input SI.HeatFlowRate Qdot_reac[n];
  input Boolean bool_eps[n];
  SI.Temperature T[n]=propsLiq.T;
  input ThermalSeparation.Units.MolarEnthalpy delta_hv[n,nSV];
  input SI.HeatFlowRate Qdot_wall[n] "heat flow rate to wall";
  input SI.MolarFlowRate Ndot_v[n] "total molar flow rate vapour";
  input SI.MolarFlowRate Ndot_v_in "total molar flow rate vapour";
  input SI.MolarFlowRate Ndot_l[n] "total molar flow rate liquid";
  input SI.MolarFlowRate Ndot_l_in "total molar flow rate vapour";
  input SI.VolumeFraction eps_liq[n]( stateSelect=StateSelect.default)
    "liquid volume fraction";
  input SI.VolumeFraction eps_vap[n](start=fill(0.99,n))
    "vapour volume fraction";

/*** entrainment ***/
  input SI.VolumeFlowRate Vdot_le[n] "liquid volume flow entrained by vapour"; // has to be supplied in extending class

/*** StartUp ***/
  input Real Ndot_source_startUp[n]
    "dummy molar flow rate to account for discharge of inert gas during startUp";
  input Boolean before_transition[n];
  parameter Boolean StartUp_CCS=falseannotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  parameter SI.Time delay_startUp=200 annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  input Real k;
  parameter Boolean smooth_startUp=false annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  parameter Boolean lowBoilingPoint[nSV]=fill(false,nSV) annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

  parameter SI.Area A "cross sectional area of the column";
  parameter SI.Length H "height of this section of the column";
  parameter Real eps "void fraction in the column";
  parameter SI.Density rho_solid[n];
  parameter SI.SpecificHeatCapacity c_solid;

/*** reaction ***/
  replaceable model Reaction = ThermalSeparation.Reaction.NoReaction constrainedby
    ThermalSeparation.Reaction.BaseReaction "model for chemical reaction"                                                                            annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

/*** initialization option ***/
   replaceable model InitOption =
      ThermalSeparation.Components.Columns.BaseClasses.Initialization.Init_T_xv_p_Ndot0
                                                                                                   constrainedby
    ThermalSeparation.Components.Columns.BaseClasses.Initialization.BaseInit
        annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

/*** thermodynamic equilibrium ***/
 replaceable model ThermoEquilibrium =
      ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid                                  constrainedby
    ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium
    "model for phase equilibrium"                                                         annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
 input Real gamma[n,nSL];

/*** StartUp ***/
  parameter Boolean considerStartUp = false
    "true if StartUp is to be considered" annotation(Dialog(tab="StartUp"));
  input SI.Pressure p_hyd[n+1] "hydraulic pressure";
  input Real omega[n];
  input Boolean startUp[n](start=fill(true,n),fixed=false);

/*** steady state ***/
  final Boolean stat "true for steady state balancing equations";

equation
  for j in 1:n loop
    c_l_star[j,:] =c_l[j,:];
    c_v_star[j,:] =c_v[j,:];
    p_sat[j,:]=p_sat_bulk[j,:];
  end for;

end BaseBalanceEquationsEq;
