within ThermalSeparation.Components.Columns;
model TrayColumn
 extends ThermalSeparation.Icons.Color.TrayColumn;

  extends ThermalSeparation.Components.Columns.BaseClasses.FeedColumn(
  A=geometry.A,
  H=geometry.H,
  eps = geometry.eps,
  rho_solid = geometry.rho_solid,
  c_solid = geometry.c_solid,
  final n= n_trays);

/*** balance equations ***/
  replaceable model BalanceEquations =
      ThermalSeparation.BalanceEquations.TrayColumn.NonEquilibrium.TwoPhaseVarState        constrainedby ThermalSeparation.BalanceEquations.TrayColumn.BaseTray
    "selection balance equations, film model and states"  annotation(choicesAllMatching=true);
  BalanceEquations balanceEquations(
                                    redeclare replaceable model HomotopyMethod =
        HomotopyMethod,
                                    redeclare replaceable package MediumLiquid =
        MediumLiquid,
                      redeclare replaceable model Reaction = Reaction,
                                                                        redeclare replaceable record Geometry =  Geometry,
                                                                        redeclare replaceable model ThermoEquilibrium =
        ThermoEquilibrium,
                                                                         redeclare replaceable model InitOption =   InitOption,
                                                                                               redeclare replaceable package MediumVapour =
        MediumVapour,
                                    final n=n,
                                              final nS=nS,
                                                                  final A=A,
                                                                                                final H=H,
                                                                                                        final eps=eps,
                                                                                                        final rho_solid=rho_solid,
                                                                                                        final c_solid=c_solid,
                                                                                                        inertLiquid=inertLiquid,
                                                                                                        inertVapour=inertVapour,
                                                                                                        propsVap=propsVap,
                                                                                                        propsVapIn=propsVapIn,
                                                                                                        stateVap=mediumVapour.state,
                                                                                                        u_v=u_v,
                                                                                                        propsLiq=propsLiq,
                                                                                                        propsLiqIn=propsLiqIn,
                                                                                                        stateLiq=mediumLiquid.state,
                                                                                                        u_l=u_l,
                                                                                                        c_v=c_v,
                                                                                                        x_v_star=x_v_star,
                                                                                                        Vdot_v_in=Vdot_v_in,
                                                                                                        Vdot_v=Vdot_v,
                                                                                                        x_upStreamIn_act=x_upStreamIn_act,
                                                                                                        x_upStreamOut_act=x_upStreamOut_act,
                                                                                                        h_upStreamIn_act=h_upStreamIn_act,
                                                                                                        h_upStreamOut_act=h_upStreamOut_act,
                                                                                                        c_l=c_l,
                                                                                                        x_l_star=x_l_star,
                                                                                                        Vdot_l_in=Vdot_l_in,
                                                                                                        Vdot_l=Vdot_l,
                                                                                                        x_downStreamIn_act=x_downStreamIn_act,
                                                                                                        x_downStreamOut_act=x_downStreamOut_act,
                                                                                                        h_downStreamIn_act=h_downStreamIn_act,
                                                                                                        h_downStreamOut_act=h_downStreamOut_act,
                                                                                                        Vdot_v_feed=Vdot_v_feed,
                                                                                                        c_v_feed=c_v_feed,
                                                                                                        h_v_feed=h_v_feed,
                                                                                                        Vdot_l_feed=Vdot_l_feed,
                                                                                                        c_l_feed=c_l_feed,
                                                                                                        h_l_feed=h_l_feed,
                                                                                                        rho_l_feed=rho_l_feed,
                                                                                                        rho_v_feed=rho_v_feed,
                                                                                                        MM_l_feed=MM_l_feed,
                                                                                                        MM_v_feed=MM_v_feed,
                                                                                                        eps_liq=eps_liq,
                                                                                                        eps_vap=eps_vap,
                                                                                                        Ndot_l_transfer=Ndot_l_transfer,
                                                                                                        Ndot_reac=Ndot_reac,
                                                                                                        Qdot_reac=Qdot_reac,
                                                                                                        bool_eps=bool_eps,
                                                                                                        delta_hv=delta_hv,
                                                                                                        Qdot_wall=Qdot_wall,
                                                                                                        Ndot_v=Ndot_v,
                                                                                                        Ndot_l=Ndot_l,
                                                                                                        Ndot_v_in=Ndot_v_in,
                                                                                                        Ndot_l_in=Ndot_l_in,
                                                                                                        Vdot_le=Vdot_le,
                                                                                                        Ndot_source_startUp=Ndot_source_startUp,
                                                                                                        mapping=mapping,
                                                                                                        gamma=activityCoeff[:].gamma,
                                                                                                        considerStartUp=considerStartUp,
                                                                                                        p_hyd=p_hyd,
                                                                                                        omega=omega,
                                                                                                        startUp=startUp,
                                                                                                        T_ref=T_ref,
                                                                                                        x_vap_liq=x_vap_liq,
                                                                                                        p_v=p_v,
                                                                                                        F=holdup.F,
                                                                                                        F_max=holdup.F_max,
                                                                                                        eps_liq_2ph=holdup.eps_liq_2ph,
                                                                                                        h=h,
                                                                                                        p_sat_bulk=p_sat_bulk,
                                                                                                        h_start=h_start,
                                                                                                        k=k,
                                                                                                       smooth_startUp=smooth_startUp,
                                                                                                        before_transition=before_transition,
                                                                                                        StartUp_CCS=StartUp_CCS,
                                                                                                        delay_startUp);

   replaceable model InitOption =
      ThermalSeparation.Components.Columns.BaseClasses.Initialization.Init_T_xv_p_Ndot0            constrainedby ThermalSeparation.Components.Columns.BaseClasses.Initialization.BaseInit
        annotation(Dialog(tab="Initialization"),choicesAllMatching=true);

 InitOption initOption(
 redeclare package MediumLiquid = MediumLiquid,
 redeclare package MediumVapour = MediumVapour,
 propsLiq = mediumLiquid.properties,   considerStartUp =   considerStartUp,
 Edot_l_transfer=Edot_l_transfer,c_l=c_l,nS=nS, mapping=mapping, n=n, nSL=nSL, inertLiquid=inertLiquid, inertVapour=inertVapour, nSV=nSV, p_v=p_v[1:n], p_v_start=p_v_start, Ndot_v_transfer=Ndot_v_transfer, x_l=x_l, x_v=x_v, x_l_start=x_l_start, x_v_start=x_v_start, n_mol_L=n_mol_L, n_mol_V=n_mol_V, x_l_star=x_l_star, x_v_star=x_v_star, x_total_start=x_total_start, T_v=T_v, T_l=T_l, T_v_start=T_v_start, T_l_start=T_l_start,rho_l=rho_l,rho_v=rho_v);
// c_A in mol_A / m3_gas ODER mol_A / m3_flssig
// rho_v in mol_ges(vap)/m3_gas oder kg_ges(vap)/m3_gas
// eps_liq in m3 Flssigkeit / m3 freies Volumen
// eps in m3 Solid / m3 gesamt

  replaceable model Reaction = ThermalSeparation.Reaction.NoReaction constrainedby ThermalSeparation.Reaction.BaseReaction
                                            "model for chemical reaction"                                                                            annotation(choicesAllMatching=true);
 Reaction reaction[n](redeclare record Geometry =  Geometry,propsLiq=mediumLiquid.properties,
       each final n= n,  c=c_l, V=A*H/n*(eps_liq), Ndot_l_transfer=Ndot_l_transfer,  gamma=activityCoeff[:].gamma,
    redeclare package MediumLiquid =   MediumLiquid);

    final parameter Integer aux[  :] = {1,3,6,10,15, 21, 28, 36, 45};

  parameter Integer n_trays(min=1)=2 "number of trays in the section";

//Initialization
  parameter SI.Height h_start[n](each max = geometry.TS) = 0.009*ones(n)
    "start value for height of the 2ph regime on the tray"                                                                      annotation(Dialog(tab="Initialization", group="Initial liquid content"));

  SI.Height h[n] "height of the 2ph region on the tray";

  replaceable record Geometry =
      ThermalSeparation.Geometry.PlateColumn.Geometry                           constrainedby ThermalSeparation.Geometry.PlateColumn.Geometry
                                                                                                                                annotation(choicesAllMatching);

  replaceable model PressureLoss =
      ThermalSeparation.PressureLoss.TrayColumn.DryHydrostatic (
      redeclare record Geometry =    Geometry)                         constrainedby ThermalSeparation.PressureLoss.TrayColumn.BasicPressureLossPlate
                                                     annotation(choicesAllMatching);

replaceable model HeatTransferWall =
      ThermalSeparation.Wall.ConstAlpha                                            constrainedby ThermalSeparation.Wall.BaseWall
                                                          annotation(choicesAllMatching=true);
    //instances of connector for heat loss
  ThermalSeparation.Interfaces.HeatPort heatPort(Qdot=-sum(heatTransferWall.Qdot_out))
    annotation (Placement(transformation(extent={{74,-10},{94,10}},
                                                                  rotation=0),
        iconTransformation(extent={{74,-10},{94,10}})));

   PressureLoss pressureLoss(redeclare model HomotopyMethod =
          HomotopyMethod,
     redeclare package MediumLiquid = MediumLiquid,
  redeclare package MediumVapour=MediumVapour,
  Vdot_l=Vdot_l, propsLiq = mediumLiquid.properties, propsVap = mediumVapour.properties,
     redeclare record Geometry = Geometry,
     final n=n,p= p_hyd, h=h, Vdot_in=Vdot_v_in,
     eps_liq_2ph=holdup.eps_liq_2ph, eps_liq=eps_liq, startUp=startUp);

  HeatTransferWall heatTransferWall(n=n,T=T, T_amb=heatPort.T, T_start = T_l_start[1], redeclare record Geometry =
                                                                                                                Geometry);
  Geometry geometry(n=n);

  replaceable model Holdup = ThermalSeparation.Holdup.TrayColumn.Stichlmair constrainedby ThermalSeparation.Holdup.TrayColumn.BaseHoldup annotation(choicesAllMatching=true);

  Holdup holdup(
    n=n, considerStartUp=considerStartUp, rho_l=rho_l, rho_v=rho_v, sigma=mediumLiquid.sigma,
    omega=omega, h=h, Vdot_v=Vdot_v,
    redeclare record Geometry =      Geometry);

   /***entrainment ***/
   Real e[n] "weight entrainment ratio";
   final parameter Real e1 = 2e-5 "value for interpolation";
   final parameter Real e2 = 1e-3 "value for interpolation";
public
   parameter Boolean entrainment = false
    "true if entrainment is to be considered";

initial equation

//    if init==Enumerations.InitializationOption.init_eq and not considerStartUp then
//      filmModel.Qdot_l_transfer=zeros(n);
//    end if;

equation
for j in 1:n-1 loop
  if startUp[j+1] then
    /*** during the startUp this requirement can be neglected; the startUp of the stage above has to be terminated before these two asserts are checked ***/
  assert(h[j]/holdup.eps_liq_2ph[j] < geometry.TS, "height of two-phase regime on tray gets larger than tray spacing!");
  assert(holdup.F_max[j]*1.2>holdup.F[j], "vapour load is higher than maximum vapour load!");
  end if;
end for;

  eps_liq = h ./ (geometry.H/n);

  /*** entrainment ***/
  for j in 1:n loop
    //Giles and Dryden: Hydrodynamics of slat trays
    if entrainment then
   Vdot_le[j] =  if eps_liq[j]< e1 then 0 else if eps_liq[j]> e2 then e[j]*Vdot_v[j]*rho_v[j] /rho_l[j] else e[j]*Vdot_v[j]*rho_v[j] /rho_l[j] *(1/(e2 - e1)*eps_liq[j] -e1/(e2-e1));
      e[j] =1e-6* 3.69e-9 * (0.073/mediumLiquid[j].sigma)*(Vdot_v[j]/A^2*rho_v[j]/1.356e-3/((geometry.TS-h[j])/0.0254))^3.2;
    else
   Vdot_le[j] = 0;
     e[j] = 0;
    end if;
  //Hunt et al: Capacity Factors in the Performance of Perforated-plate Columns: irgendwie ergibt das sehr groe Werte fr Vdot_le!
  //  Vdot_le[j] =  if eps_liq[j]<2e-5 then 0 else e[j]*Vdot_v[j]*MM_v[j]/MM_l[j]*eps_liq[j]^(1/50);
 // Vdot_le[j] =  if eps_liq[j]<2e-5 then 0 else e[j]*Vdot_v[j]*MM_v[j]/MM_l[j];
  // e[j] = 0.22 * (0.073/mediumLiquid[j].sigma)*(Vdot_v[j]/A/0.3048/((geometry.TS-geometry.TS*eps_liq[j])/0.0254))^3.2;
    end for;

      /*** variables inherited from the base class: ***/
      Vdot_v = pressureLoss.Vdot;
         p_v_in=pressureLoss.p_v_in;
        Qdot_wall = heatTransferWall.Qdot;
   Vdot_l = holdup.Vdot;
   Ndot_reac = reaction.Ndot;
   Qdot_reac = reaction.deltaH_R;

/*** variables determined in the film model: ***/
     x_v=balanceEquations.x_v;
     Ndot_v_transfer=balanceEquations.Ndot_v_transfer;
     Edot_l_transfer= balanceEquations.Edot_l_transfer;
     Edot_v_transfer= balanceEquations.Edot_v_transfer;
     T_star = balanceEquations.T_star;

                                      annotation (extent=[-70,-100; -50,-80], Diagram(graphics),
               Documentation(revisions="<html>
<ul>
 
<table border>
 
<tr>
  <th> created by </th>
  <td> <a href=\"mailto:karin.dietl@tu-harburg.de\">Karin Dietl</a> &amp; <a href=\"mailto:andreas.joos@tu-harburg.de\">Andreas Joos</a></td>
</tr>
<tr>
  <th> creation date </th>
  <td> 01.01.2009 </td>
</tr>
 
<tr>
  <th> revised by </th>
  <td> nobody so far</td>
</tr>
<tr>
  <th> last revision </th>
  <td> this is an alpha version... </td>
</tr>
<tr>
  <th> based on </th>
  <td> </td>
</tr>
</table>
 
 
</ul> 
 
</html>", info="<html>
<p>This model discribes a tray column. </p>
<p>Base classes for the following classes are instantiated her:</p>
<p><ul>
<li>liquid holdup</li>
<li>pressure loss </li>
<li>film model (<a href=\"Modelica://ThermalSeparation.FilmModel.TrayColumn.BaseFilmTray\">BaseFilmTray</a>) </li>
<li>geometry</li>
<li>heat transfer to wall</li>
</ul></p>
<p><h4>Volume flow rate of the liquid leaving the tray </h4></p>
<p>Liquid leaves the tray if the height of the two-phase regime on the tray, h, gets higher than the height of the weir, h_w. The volume flow rate of the liquid, Vdot_l, is then proportional to the height over weir, h_ow = h-h_w and is calculated using the following formula: </p>
<p>Vdot_l = l_w &middot; eps_liq_2ph &middot; ((h / eps_liq_2ph - h_w) &middot; g1/3 / 1.45)3/2 </p>
<p>where l_w is the weir length and eps_liq_2ph the liquid fraction in the two-phase regime on the tray. This equation was obtained from Stichlmair [1]. However this equation would yield a negative volume flow rate, if h &LT; h_w. Since this is not possible, Vdot_l is set to zero if h &LT; h_w (which may be the case for the start-up of the column). It is supposed that there is only liquid leaving the tray, if the height of the two-phase regime is high enough; &QUOT; raining &QUOT; through the holes in the tray (which occurs if the vapour load is too small) is not considered. </p>
<p><h4>Operating range </h4></p>
<p>A very important point for tray columns is the operating range. There exists a minimum and a maximum vapour load as well as a minimum and maximum liquid load. In order to discribe the minimum vapour flow not the vapour flow Vdot_v itself is used but the vapour load F (which is defined as F = wV &middot; &rho;V wV is the superficial velocity) or the vapour load Fh, where the velocity is the velocity in the holes. In theory the F-values are different for all elements of the column, however the F-values are calculated only once for each section using the inlet conditions for vapour and liquid. This is acceptable since during normal operation the F-values don&apos;t vary a lot over one section and for conditions like for example start-up operation the equations to calulate maximum and minimum load are not valid anyway. Also the F-values shall only give an idea of the operation range and in any case it is suggested to stay well in the operating range. </p>
<p><h4>&nbsp; &nbsp; Minimum vapour load </h4></p>
<p>&nbsp; &nbsp; A minimum vapour load exists for sieve trays, since if the vapour load is too small, liquid is going to rain through the plates. To avoid this Ruff [2] found out that the vapour load must be at least</p>
<p>&nbsp; &nbsp; Fh,min,Ruff &ge; (0.37 &middot; dh &middot; (&rho;L - &rho;V)5/4 / &rho;V1/4)0.5 </p>
<p>&nbsp; &nbsp; where dh is the diameter of the holes in the tray and &rho;L and &rho;V are the densities of the liquid and the vapour respectively. However if the holes in the trays are rather small (around 2 mm - 3 mm) raining is not the major problem, but the fact that vapour is only passing a part of the tray. This can be avoided the following holdes: </p>
<p>&nbsp; &nbsp; Fh,min,Mersmann &ge; (2 &middot; &AMP;sigma / dh)0.5 </p>
<p>&nbsp; &nbsp; This formula was found by Mersmann [3]. So the minimum vapour load was defined to be Fh,min = min(Fh,min,Ruff ,Fh,min,Mersmann). The model allows for vapour loads which are smaller than the minimum vapour load (so no <code>assert</code> is used in the code), however one has to bear in mind that raining and maldistribution of the gas are not modelled and the results for Fh &LT; Fh,min are not very reliable. </p>
<p><h4>&nbsp; &nbsp; Maximum vapour load </h4></p>
<p>&nbsp; &nbsp; If the vapour load is too high, <b>entrainment</b> occurs, i.e. vapour blows the liquid out of the column. An equation from Stichlmair [1] was taken in order to calulate the maximum vapour load: </p>
<p>&nbsp; &nbsp; Fmax = 2.5 &middot; (&phi;2 &middot; &sigma; &middot; (&rho;L - &rho;V) &middot; g)1/4 &middot; (100 Vdot_l/Vdot_v)0.06 / (1 - 10 Vdot_l/Vdot_v)0.5 </p>
<p>&nbsp; &nbsp; Again the compliance with this rule is not ensured via an <code>assert</code>, but the user is advised to check that F &LT; Fmax. If this is not the case the equation used to calculate Vdot_l is not longer valid. </p>
<p><h4>&nbsp; &nbsp; Minimum and maximum liquid load </h4></p>
<p>&nbsp; &nbsp; Theoretically the liquid load can get very small, however this will not be very efficient, so it is recommanded to have a minimum height over weir of 5 mm. The liquid flows through the column due to the earth gravity. Therefore the capacity is limited. However up to now no correlations are implemented. </p>
<p><h4>Literature </h4></p>
<p>[1] Stichlmair: Dimensionierung des Gas/Fl&uuml;ssigkeits-Kontaktapparates Bodenkolonne, Chem.-Ing.-Tech. 50 (1978), Nr. 4, p. 281-284 </p>
<p>[2] Ruff et al., Chem.-Ing.Tech. 48 (1976) Nr. 9, p. 759-764 </p>
<p>[3] Mersmann, Chem.-Ing.Tech. 35 (1963) Nr. 2, p. 103-107 </p>
</html>"),
    Placement(transformation(extent={{-70,-100},{-50,-80}}, rotation=0)),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}), graphics));
end TrayColumn;
