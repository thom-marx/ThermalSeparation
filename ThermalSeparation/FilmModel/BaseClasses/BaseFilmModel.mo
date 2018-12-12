within ThermalSeparation.FilmModel.BaseClasses;
partial model BaseFilmModel
  "base class for film model where no source term for the reaction exists"
  replaceable model Reaction = ThermalSeparation.Reaction.NoReaction                               constrainedby
    ThermalSeparation.Reaction.BaseReaction "model for chemical reaction" annotation(choicesAllMatching=true,Dialog(
      tab="Propagated from Column",
      group=
          "These variables are propagated from the column model and do not have to be set by the user!",
      enable=false));

  replaceable package MediumLiquid = Media.BaseMediumLiquid annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  replaceable package MediumVapour = Media.BaseMediumVapour annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  parameter Boolean inertVapour[nSV] = fill(false,nSV)
   annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  parameter Boolean inertLiquid[nSL] = fill(false,nSL) annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  final parameter Integer nSL = MediumLiquid.nSubstance;
  final parameter Integer nSV = MediumVapour.nSubstance;
  parameter Integer n(min=1) annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  parameter Integer nS(min=2)
   annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  parameter Integer mapping[nS,2]= {{i,i} for i in 1:nS} annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  parameter Boolean enableDialog=true annotation(Dialog(enable=false));

    /*** vapour properties ***/
    input MediumVapour.ThermodynamicProperties[n] propsVap annotation(HideResult=true);
input SI.MoleFraction x_v_star[ n,nSV];
  input MediumVapour.ThermodynamicState[n] stateVap  annotation(HideResult=true);
  input SI.MoleFraction x_v_in[nSV]  annotation(HideResult=true);
  input SI.Concentration c_v_star[n,nSV];
  input SI.MoleFraction x_l[n, nSL]  annotation(HideResult=true);
output SI.MoleFraction x_v[n, nSV];
  input SI.Pressure p_sat[n,nSL]  annotation(HideResult=true);
  /*** liquid properties ***/
  input MediumLiquid.ThermodynamicProperties[n] propsLiq  annotation(HideResult=true);
  input MediumLiquid.ThermodynamicState[n] stateLiq  annotation(HideResult=true);
input SI.MoleFraction x_l_star[ n,nSL]  annotation(HideResult=true);
  input SI.Concentration c_l[ n,nSL]  annotation(HideResult=true);
  input SI.Concentration c_l_star[ n,nSL]  annotation(HideResult=true);
  input Real eta_comp[n,nSL]  annotation(HideResult=true);
  input Real gamma[n,nSL];

  input SI.MoleFraction x_vap_liq[n,nS];

//Variables upStream
  input SI.Concentration c_v_in[ nSV]  annotation(HideResult=true);
  input SI.Concentration c_v[ n,nSV];
  input SI.VolumeFlowRate Vdot_v_in( nominal=1e-2)  annotation(HideResult=true);
  input SI.VolumeFlowRate Vdot_v[ n](nominal=fill(1e-2,n))  annotation(HideResult=true);
  input SI.VolumeFlowRate Vdot_l[n]  annotation(HideResult=true);
  input Real eps_liq[n]  annotation(HideResult=true);

  input SI.MolarFlowRate Ndot_l_transfer[          n,nSL](start=fill(0.1,n,nSL))  annotation(HideResult=true);
output SI.MolarFlowRate Ndot_v_transfer[          n,nSV](start=fill(-0.1,n,nSV))
    "positive when entering the bulk phase";

     /*** energy balance ***/
output SI.HeatFlowRate Edot_l_transfer[     n];
output SI.HeatFlowRate Edot_v_transfer[  n];

  parameter SI.Temperature T_ref   annotation(HideResult=true,Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  input SI.Pressure p_hyd[n+1]  annotation(HideResult=true);
  input SI.Pressure p_v[n+1]  annotation(HideResult=true);

output SI.Temperature T_star[n](start=fill(273.15,n));

   /*** StartUp ***/
  parameter Boolean considerStartUp = false
  annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  input Real omega[n];
  input Boolean startUp[n];
  input Boolean before_transition[n];
  parameter Boolean StartUp_CCS=false;
  parameter SI.Time delay_startUp=200 annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  input Real k;
  parameter Boolean smooth_startUp=false annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));
  parameter Boolean lowBoilingPoint[nSV]=fill(false,nSV) annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

  replaceable record BaseGeometry =
      ThermalSeparation.Geometry.BasicGeoPackedColumn
  constrainedby ThermalSeparation.Geometry.BasicGeometry;
//has to be provided by the extending class: true if the extending class is an equilibrium model, false otherwise
parameter Boolean EQ;

Real rel_dev_liq[n,nSL] "relative deviation from equilibrium, liquid";
Real rel_dev_vap[n,nSV] "relative deviation from equilibrium, vapour";
Real max_rel_dev_liq = max(rel_dev_liq)
    "maximum relative deviation from equilibrium, liquid";
Real max_rel_dev_vap = max(rel_dev_vap)
    "maximum relative deviation from equilibrium, vapour";

  /*** Homotopy ***/
  // Boolean useHomotopy;
  replaceable model HomotopyMethod =
      ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.NoHomotopy
                       constrainedby
    ThermalSeparation.Components.Columns.BaseClasses.Initialization.Homotopy.BaseHomotopy;

HomotopyMethod homotopyMethod(nS=nS,n=n,nSV=nSV,nSL=nSL);

protected
  SI.Temperature T_l[n] = propsLiq.T;
  SI.Temperature T_v[n] = propsVap.T;
  SI.MolarInternalEnergy u_v[n] = propsVap.u;

  SI.MolarFlowRate Ndot_v_interface[n,nSV](start=fill(-0.1,n,nSV))
    "positive when going from interface to vapour film";
 SI.MolarFlowRate Ndot_l_interface[n,nSL](start=fill(0.1,n,nSL))
    "positive when going from interface to liquid film";
  SI.HeatFlowRate Edot_v_interface[n];
  SI.HeatFlowRate Edot_l_interface[ n];

equation
  for j in 1:n loop
//  x_l[j,:] = propsLiq[j].x;
//  x_v[j,:] = propsVap[j].x;
  end for;

/*** molar flux over phase boundary zero for inert substances ***/
  for j in 1:n loop
    for i in 1:nSV loop
      if inertVapour[i] then
      Ndot_v_interface[j,i] = 0;
      end if;
    end for;
    for i in 1:nSL loop
      if inertLiquid[i] then
      Ndot_l_interface[j,i]= 0;
      end if;
    end for;
  end for;

/*** balance equations at phase boundary: always steady-state ***/
  - Edot_l_interface[:] - Edot_v_interface[:] = zeros(n);

  for i in 1:nS loop
    Ndot_v_interface[:,mapping[i,1]] + Ndot_l_interface[:,mapping[i,2]]  = zeros(n);
  end for;

   for j in 1:n loop
  for i in 1:nSL loop
    rel_dev_liq[j,i] = (abs(-x_l[j,i] + x_l_star[j,i]))/max(1e-9,x_l[j,i]);
  end for;
  for i in 1:nSV loop
    rel_dev_vap[j,i] = (abs(-x_v[j,i] + x_v_star[j,i]))/max(1e-9,x_v[j,i]);
  end for;
  end for;

  annotation (Documentation(info="<html>
<p>Base class for a film model. It provides the interface for all film models.</p>
<p>Additionally the balance equations at the phase boundary are written here: since the phase boundary is always calculated steady-state, the molar flow rates and the energy flow rates sum up to zero. Also the equations for the phase boundary molar flow rates of inert components are written in this class, namely: </p>
<p>N_flow_inert_PB = 0</p>
</html>"));
end BaseFilmModel;
