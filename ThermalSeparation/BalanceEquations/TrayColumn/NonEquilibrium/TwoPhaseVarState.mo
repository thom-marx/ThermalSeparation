within ThermalSeparation.BalanceEquations.TrayColumn.NonEquilibrium;
model TwoPhaseVarState "phases balanced seperately, states optional"

extends
    ThermalSeparation.BalanceEquations.Base.NonEquilibrium.BaseTwoPhaseVarState;
extends ThermalSeparation.BalanceEquations.TrayColumn.BaseTray;

parameter Boolean EQ=filmModel.EQ;

  input Real eps_liq_2ph[n]
    "liquid fraction in the two-phase area on the plate";
  input ThermalSeparation.Units.F_Factor F[n]
    "F-factor based on the active area";
  input ThermalSeparation.Units.F_Factor F_max[n] "maximum vapour load";
  input SI.Height h[n];

/*** geometry ***/
  replaceable record Geometry =
      ThermalSeparation.Geometry.PlateColumn.Geometry                constrainedby ThermalSeparation.Geometry.PlateColumn.Geometry
                                                    "column geometry"          annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

 /*** film model ***/
   replaceable model FilmModel =
       ThermalSeparation.FilmModel.TrayColumn.MS (redeclare replaceable model
        StateSelection =
      ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.StateSelection2)       constrainedby ThermalSeparation.FilmModel.TrayColumn.BaseFilmTray(
     redeclare replaceable package MediumLiquid =  MediumLiquid,
         redeclare replaceable model Reaction =  Reaction,
         redeclare replaceable package MediumVapour =  MediumVapour)
    "heat and mass transfer mechanism across phase boundary and state selection"  annotation(choicesAllMatching=true);
 FilmModel filmModel(redeclare replaceable model Reaction = Reaction,
 propsLiq = propsLiq, propsVap = propsVap,
   redeclare record Geometry =  Geometry,p_hyd=p_hyd, p_v=p_v, T_ref=T_ref,nS=nS, mapping=mapping,
   redeclare package MediumLiquid =   MediumLiquid,x_v_star=x_v_star, x_l_star=x_l_star,
   redeclare package MediumVapour =   MediumVapour,                                                                                                    inertVapour=inertVapour, inertLiquid=inertLiquid, final n=n,
   Vdot_l=Vdot_l,  c_l=c_l, c_l_star=c_l_star,  c_v_in=c_v_in, c_v=c_v, Vdot_v_in=Vdot_v_in, Vdot_v=Vdot_v, considerStartUp=considerStartUp,  omega=omega, startUp=startUp,
  stateLiq=stateLiq, stateVap=stateVap, Ndot_l_transfer=Ndot_l_transfer,
  eps_liq=eps_liq, eta_comp=propsLiq.eta_comp, x_l=x_l, p_sat=p_sat,x_v_in=x_v_in,gamma=gamma,
  redeclare final model ThermoEquilibrium =     ThermoEquilibrium,
  c_v_star=c_v_star,x_vap_liq=x_vap_liq,F=F, F_max=F_max, eps_liq_2ph=eps_liq_2ph,h=h,redeclare replaceable model
                        HomotopyMethod =
        HomotopyMethod,
          k=k,
          smooth_startUp=smooth_startUp,
          before_transition=before_transition,
          StartUp_CCS=StartUp_CCS,
          delay_startUp);

equation
     x_v=filmModel.x_v;
     Ndot_v_transfer=filmModel.Ndot_v_transfer;
     Edot_l_transfer= filmModel.Edot_l_transfer;
     Edot_v_transfer= filmModel.Edot_v_transfer;
     T_star = filmModel.T_star;

initial equation
  h = h_start;

end TwoPhaseVarState;
