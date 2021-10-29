within ThermalSeparation.Components.ColumnsNoIndex;
model StructuredPackedColumn_external
     extends ThermalSeparation.Icons.Color.PackedColumn;
  extends ThermalSeparation.Components.ColumnsNoIndex.BaseClasses.FeedColumn(
    A=geometry.A,
    H=geometry.H,
    eps=geometry.eps,
    rho_solid=geometry.rho_solid,
    c_solid=geometry.c_solid,
    final n=n_elements,
    final EQ=filmModel.EQ);

   replaceable model InitOption =
      ThermalSeparation.Components.ColumnsNoIndex.BaseClasses.Initialization.Init_T_xv_p_Ndot0
    constrainedby ThermalSeparation.Components.ColumnsNoIndex.BaseClasses.Initialization.BaseInit
        annotation(Dialog(tab="Initialization"),choicesAllMatching=true);

 InitOption initOption(
 redeclare package MediumLiquid = MediumLiquid,
 redeclare package MediumVapour = MediumVapour,
 propsLiq = mediumLiquid.properties,   considerStartUp =   considerStartUp,
 Edot_l_transfer=Edot_l_transfer,c_l=c_l,nS=nS, mapping=mapping, n=n, nSL=nSL, inertLiquid=inertLiquid, inertVapour=inertVapour, nSV=nSV, p_v=p_v[1:n], p_v_start=p_v_start, Ndot_v_transfer=Ndot_v_transfer, x_l=x_l, x_v=x_v, x_l_start=x_l_start, x_v_start=x_v_start, n_mol_L=n_mol_L, n_mol_V=n_mol_V, x_l_star=x_l_star, x_v_star=x_v_star, x_total_start=x_total_start, T_v=T_v, T_l=T_l, T_v_start=T_v_start, T_l_start=T_l_start);
// c_A in mol_A / m3_gas ODER mol_A / m3_flüssig
// rho_v in mol_ges(vap)/m3_gas oder kg_ges(vap)/m3_gas
// eps_liq in m3 Flüssigkeit / m3 freies Volumen
// eps in m3 Solid / m3 gesamt

  replaceable model Reaction = ThermalSeparation.Reaction.NoReaction constrainedby ThermalSeparation.Reaction.BaseReaction
                                            "model for chemical reaction"                                                                            annotation(choicesAllMatching=true);
 Reaction reaction[n](redeclare record Geometry =  Geometry,propsLiq=mediumLiquid.properties,
       each final n= n,  c=c_l, V=A*H/n*(eps_liq), Ndot_l_transfer=Ndot_l_transfer,  gamma=activityCoeff[:].gamma,
    redeclare package MediumLiquid =   MediumLiquid);

  replaceable model FilmModel =
      ThermalSeparation.FilmModel.StructuredPackedColumn.MS (
    redeclare replaceable package MediumLiquid =  MediumLiquid,
        redeclare replaceable model Reaction =  Reaction,
        redeclare replaceable package MediumVapour =  MediumVapour)                 constrainedby ThermalSeparation.FilmModel.StructuredPackedColumn.BaseFilmPacked(
    redeclare replaceable package MediumLiquid =  MediumLiquid,
        redeclare replaceable model Reaction =  Reaction,
        redeclare replaceable package MediumVapour =  MediumVapour)
    "heat and mass transfer mechanism across phase boundary"                                                                                                     annotation(choicesAllMatching=true);
FilmModel filmModel(redeclare replaceable model Reaction = Reaction,
propsLiq = mediumLiquid.properties, propsVap = mediumVapour.properties,
  redeclare record Geometry =  Geometry,p_hyd=p_hyd, p_v=p_v, T_ref=T_ref,nS=nS, mapping=mapping,
  redeclare package MediumLiquid =   MediumLiquid,x_v_star=x_v_star, x_l_star=x_l_star,
  redeclare package MediumVapour =   MediumVapour,                                                                                                    inertVapour=inertVapour, inertLiquid=inertLiquid, final n=n,
  Vdot_l=Vdot_l,  c_l=c_l, c_l_star=c_l_star,  c_v_in=c_v_in, c_v=c_v, Vdot_v_in=Vdot_v_in, Vdot_v=Vdot_v, considerStartUp=considerStartUp,  omega=omega, startUp=startUp,
 stateLiq=mediumLiquid.state, stateVap=mediumVapour.state, Ndot_l_transfer=Ndot_l_transfer,
 eps_liq=eps_liq, eta_comp=mediumLiquid.eta_comp, x_l=x_l, p_sat=p_sat,x_v_in=x_v_in,gamma=activityCoeff[:].gamma,
 redeclare final model ThermoEquilibrium =     ThermoEquilibrium,
 c_v_star=c_v_star,x_vap_liq=x_vap_liq);

  final parameter Integer aux[  :] = {1,3,6,10,15, 21, 28, 36, 45};

  parameter Integer n_elements(min=1)=2
    "number of discrete elements or number of equilibrium stages (for equilibrium film model)";

//Initialization
  parameter Boolean wettedInitial = true
    "true if package material is already wetted (i.e. if liquid enters the column there is immediately a liquid stream leaving the column)"
                                                                                                        annotation(Dialog(tab="Initialization", group="Initial liquid content"));
   parameter Real eps_liq_start = 0.06
    "start value for liquid content if it is not exactly wetted but with more or less liquid"  annotation(Dialog(enable=not wettedInitial, tab="Initialization", group="Initial liquid content"));

  replaceable record Geometry =
      ThermalSeparation.Geometry.StructuredPackedColumn.Geometry                constrainedby ThermalSeparation.Geometry.StructuredPackedColumn.Geometry
    "column geometry"                                                                           annotation(choicesAllMatching);

   replaceable model PressureLoss =
      ThermalSeparation.PressureLoss.StructuredPackedColumn.Stichlmair ( redeclare record
               Geometry =                                                                             Geometry)                         constrainedby ThermalSeparation.PressureLoss.StructuredPackedColumn.BasicPressureLossPacked
    "pressure loss correlation"                                                  annotation(choicesAllMatching);

  PressureLoss pressureLoss(
    redeclare package MediumLiquid = MediumLiquid,
  redeclare package MediumVapour=MediumVapour,
  Vdot_l=Vdot_l, propsLiq = mediumLiquid.properties, propsVap = mediumVapour.properties,
    redeclare record Geometry = Geometry, final n=n,p= p_hyd,
    Vdot_in=Vdot_v_in, hu_dyn = holdup.hu_dyn, eps_liq=eps_liq);

  replaceable model HeatTransferWall =
      ThermalSeparation.Wall.ConstAlpha                                            constrainedby ThermalSeparation.Wall.BaseWall
    "heat transfer mechanism between bulk and wall"                                                                                                annotation(choicesAllMatching=true);
    //instances of connector for heat loss
  ThermalSeparation.Interfaces.HeatPort heatPort(Qdot=-sum(heatTransferWall.Qdot_out))
    annotation (Placement(transformation(extent={{74,-10},{94,10}},
                                                                  rotation=0),
        iconTransformation(extent={{74,-10},{94,10}})));

  replaceable model Holdup =
      ThermalSeparation.Holdup.StructuredPackedColumn.StichlmairStat              constrainedby ThermalSeparation.Holdup.StructuredPackedColumn.BaseHoldup
    "calculation of holdup inside the column"                                                                                                     annotation(choicesAllMatching);
  Holdup holdup(
      redeclare package MediumLiquid = MediumLiquid,
  redeclare package MediumVapour=MediumVapour,
  propsLiq = mediumLiquid.properties, propsVap = mediumVapour.properties,p_v= p_v,
    final n=n,
    redeclare record Geometry=Geometry,
    eps_liq=eps_liq);

  Geometry geometry(n=n);

  HeatTransferWall heatTransferWall(
    final n=n,T=T, T_amb=heatPort.T, T_start = T_l_start[1],
    redeclare record Geometry = Geometry);

  Real F[n]= Vdot_v/geometry.A .* rho_v.^0.5
    "F-Factor to check operating range, unit: Pa^0.5";
  Real liquidLoad[n] = Vdot_l/geometry.A*3600
    "liquid load to check operating range, unit: m3/m2/h";

initial equation
   if wettedInitial then
     eps_liq = holdup.hu_stat+fill(Modelica.Constants.eps,n);
   else
     eps_liq = fill(eps_liq_start,n);
   end if;

//    if init==Enumerations.InitializationOption.init_eq and not considerStartUp then
//      filmModel.Qdot_l_transfer=zeros(n);
//    end if;

//der(eps_liq)=zeros(n);
equation

     for j in 1:n loop
      Vdot_l[j] = max(0,holdup.Vdot[j]);
      end for;
        Vdot_v = pressureLoss.Vdot;
          p_v_in = pressureLoss.p_v_in;
        Qdot_wall = heatTransferWall.Qdot;
  Vdot_le = zeros(n);
     Ndot_reac = reaction.Ndot;
   Qdot_reac = reaction.deltaH_R;

/*** variables determined in the film model: ***/

x_v=filmModel.x_v;
Ndot_v_transfer=filmModel.Ndot_v_transfer;
Edot_l_transfer= filmModel.Edot_l_transfer;
Edot_v_transfer= filmModel.Edot_v_transfer;
T_star = filmModel.T_star;

      annotation (extent=[-70,-100; -50,-80], Diagram(graphics),
    Documentation(info="<html>
<p>Specialized class to describe a packed column. Base classes for the following classes are instantiated her:</p>
<p><ul>
<li>liquid holdup</li>
<li>pressure loss </li>
<li>film model (<a href=\"Modelica://ThermalSeparation.FilmModel.StructuredPackedColumn.BaseFilmPacked\">BaseFilmPacked</a>) </li>
<li>geometry</li>
<li>heat transfer to wall</li>
</ul></p>
<p><br/>This class of a structured packed column also supplies the initial equation for the liquid content in the column. </p>
<p>This value can be specified by the user, or the boolean parameter <i>wetted</i> is set to true, which means that the packing material at simulation start is completely wetted. For a fast simulation it is often helpful to initialize with a liquid content which is slightly above the liquid content for a wetted column.</p>
</html>"),     Documentation(revisions="<html>
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
Model eines Packungskolonnenabschnitts mit einem Dampf- und einem Fluessigstrom.
</html>"),
    Placement(transformation(extent={{-70,-100},{-50,-80}}, rotation=0)),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}),
         graphics));
end StructuredPackedColumn_external;
