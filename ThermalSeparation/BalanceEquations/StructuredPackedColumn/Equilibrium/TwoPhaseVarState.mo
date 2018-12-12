within ThermalSeparation.BalanceEquations.StructuredPackedColumn.Equilibrium;
model TwoPhaseVarState
  "Equilibrium: both phases balanced together, states optional"

extends
    ThermalSeparation.BalanceEquations.Base.Equilibrium.BaseTwoPhaseVarState;
extends
    ThermalSeparation.BalanceEquations.StructuredPackedColumn.BaseStructured;

parameter Boolean EQ=true;
Real K[n,nS];

/*** geometry ***/
  replaceable record Geometry =
      ThermalSeparation.Geometry.StructuredPackedColumn.Geometry                constrainedby
    ThermalSeparation.Geometry.StructuredPackedColumn.Geometry
    "column geometry"          annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

  /*** film model ***/
    replaceable model FilmModel =
        ThermalSeparation.FilmModel.StructuredPackedColumn.MS (redeclare
        replaceable model StateSelection =
       ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.StateSelection2)
                                                                                                        constrainedby
    ThermalSeparation.FilmModel.StructuredPackedColumn.BaseFilmPacked(
       redeclare replaceable package MediumLiquid =  MediumLiquid,
           redeclare replaceable model Reaction =  Reaction,
          redeclare replaceable package MediumVapour =  MediumVapour,
          k=k,
          smooth_startUp=smooth_startUp,
          before_transition=before_transition,
          StartUp_CCS=StartUp_CCS,
          delay_startUp)
    "heat and mass transfer mechanism across phase boundary and state selection"   annotation(choicesAllMatching=true);

  ThermoEquilibrium thermoEquilibrium[n](x_vap_liq=x_vap_liq,
    each nS=nS,  each mapping =                                    mapping,
      redeclare replaceable package MediumVapour =  MediumVapour,
      redeclare replaceable package MediumLiquid =    MediumLiquid,
       p=p_v[1:n], T=T_star, x_v=x_v_star, x_l=x_l_star, p_sat=p_sat,  v_v=MM_v./rho_v);

HomotopyMethod homotopyMethod(n=n,nS=nS,nSL=nSL,nSV=nSV);

equation
 /*** balance equations at phase boundary: always steady-state ***/
   - Edot_l_transfer[:] - Edot_v_transfer[:] = zeros(n);

   for i in 1:nS loop
     Ndot_v_transfer[:,mapping[i,1]] + Ndot_l_transfer[:,mapping[i,2]]  = zeros(n);
   end for;

       /*** energy balance at the phase boundary ***/
 for j in 1:n loop

 propsVap[j].T = T_star[j];
 propsLiq[j].T = T_star[j];
 end for;

    x_v_star = x_v;

 /*** thermodynamic equilibrium at phase boundary ***/
if homotopyMethod.bool_K then
    for j in 1:n loop
   for i in 1:nS loop
    K[j,i] = homotopy(actual=thermoEquilibrium[j].K[i],simplified=homotopyMethod.K[i]);
   end for;

end for;
else

  for j in 1:n loop
   for i in 1:nS loop
    K[j,i] = thermoEquilibrium[j].K[i];
   end for;

end for;
end if;

 for j in 1:n loop
    for i in 1:nS loop
     x_v_star[j,mapping[i,1]]= K[j,i] *x_l_star[j,mapping[i,2]];
     // K[j,i] = thermoEquilibrium[j].K[i];
    end for;

    for i in 1:nSL loop
    x_l_star[j,i] = x_l[j,i];
    end for;
 end for;

initial equation

   if wettedInitial then
     eps_liq = hu_stat+fill(Modelica.Constants.eps,n);
   else
     eps_liq = fill(eps_liq_start,n);
   end if;

end TwoPhaseVarState;
