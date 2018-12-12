within ThermalSeparation.FilmModel.BaseClasses;
model EquilibriumEffectiveDiff
  "equilibrium model using effective diffusivity method"
  extends ThermalSeparation.FilmModel.BaseClasses.BaseFilmModel(
                         final EQ=true);

Real K[n,nS] "equilibrium constant" annotation(Dialog(__Dymola_initialDialog=true));

  SI.Density rho_v[n]=propsVap.rho;
  SI.MolarMass MM_v[n] = propsVap.MM;
   replaceable model ThermoEquilibrium =
      ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid                                  constrainedby
    ThermalSeparation.PhaseEquilibrium.BasePhaseEquilibrium
       annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

  ThermoEquilibrium thermoEquilibrium[n](x_vap_liq=x_vap_liq,
    each nS=nS,  each mapping =                                    mapping,
      redeclare replaceable package MediumVapour =  MediumVapour,
      redeclare replaceable package MediumLiquid =    MediumLiquid,
       p=p_v[1:n], T=T_star, x_v=x_v_star, x_l=x_l_star, p_sat=p_sat,  v_v=MM_v./rho_v);

    /*** film is calculated steady-state, no source term due to reaction ***/
  parameter Real factorHT = 10000
    "factor for heat transfer corresponding to alpha*area";

      /*** state variables ***/
  replaceable model StateSelection =
      ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.BaseStateSelectionNoneq
                                                                                                        constrainedby
    ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.BaseStateSelectionNoneq
                                                                                annotation(choicesAllMatching,Dialog(enable=enableDialog));
  StateSelection stateSelection(
     redeclare replaceable package MediumVapour =  MediumVapour,
    redeclare replaceable package MediumLiquid =   MediumLiquid,
  n=n, propsLiq=propsLiq, propsVap=propsVap, c_l=c_l, p_v= p_v[1:n]);
  BaseGeometry baseGeometry(n=n);

  parameter ThermalSeparation.Units.CoefficentOfMassTransfer k_l = 1;
  parameter ThermalSeparation.Units.CoefficentOfMassTransfer k_v = 100;
/*** to be provided by extending class ***/
  SI.Area A_I[n] "interfacial area";
  SI.MoleFraction max_rel_error_liq = max(vector_liq);
  SI.MoleFraction max_rel_error_vap = max(vector_vap);

protected
   Real vector_liq[ n,nSL];
   Real vector_vap[ n,nSV];

       SI.HeatFlowRate Qdot_l_transfer[n];
 SI.HeatFlowRate Qdot_v_transfer[  n];

   SI.MoleFraction x_v_star_eq[n,nSV];

equation
    /*** film is calculated steady-state ***/
  if homotopyMethod.bool_Edot_inter and homotopyMethod.useHomotopy then

    Edot_v_transfer=homotopy(actual=Edot_v_interface,simplified=homotopyMethod.Edot_v_inter);
    Edot_l_transfer=homotopy(actual=Edot_l_interface,simplified=homotopyMethod.Edot_l_inter);

else

    Edot_v_transfer=Edot_v_interface;
    Edot_l_transfer=Edot_l_interface;

end if;

 if homotopyMethod.bool_Ndot_inter and homotopyMethod.useHomotopy then

  Ndot_v_transfer=homotopy(actual=Ndot_v_interface,simplified=homotopyMethod.Ndot_v_inter);
  Ndot_l_transfer=homotopy(actual=Ndot_l_interface,simplified=homotopyMethod.Ndot_l_inter);

else

  Ndot_v_transfer=Ndot_v_interface;
  Ndot_l_transfer=Ndot_l_interface;

end if;

//   Ndot_v_interface=Ndot_v_transfer;
//   Ndot_l_interface=Ndot_l_transfer;
//   Edot_v_interface=Edot_v_transfer;
//   Edot_l_interface=Edot_l_transfer;

 for j in 1:n loop
  for i in 1:nSL loop
    vector_liq[j,i] = (abs(-x_l[j,i] + x_l_star[j,i]))/max(1e-9,x_l[j,i]);
  end for;
  for i in 1:nSV loop
    vector_vap[j,i] = (abs(-x_v[j,i] + x_v_star[j,i]))/max(1e-9,x_v[j,i]);
  end for;
  end for;

 /*** thermodynamic equilbrium also in the bulk phases ***/
   for j in 1:n loop
    for i in 1:nSV loop
           if considerStartUp and startUp[j] then
         Ndot_v_transfer[j,i] = if j==1 then -Vdot_v_in*c_v_in[i] else -Vdot_v[j-1]*c_v[j-1,i];
       else
       Ndot_v_transfer[j,i] =  - A_I[j]*k_v* propsVap[j].rho/propsVap[j].MM*(x_v[j,i] - x_v_star[j,i]);

           end if;
     end for;
     //Liquid side
    for i in 1:nSL loop
          Ndot_l_transfer[j,i] = - A_I[j]*k_l*1/propsLiq[j].v*(x_l[j,i] - x_l_star[j,i]);// - geometry.H/n * factorMT_liq * (x_l[j,i] - x_l_star[j,i]);

   end for;
   end for;

      /*** energy balance at the phase boundary ***/
for j in 1:n loop

    Qdot_l_transfer[j] = baseGeometry.H/n * factorHT* (T_star[j] - T_l[j]);
    Qdot_v_transfer[j] = baseGeometry.H/n * factorHT*(T_star[j] - T_v[j]);
    Edot_l_transfer[j] =   Qdot_l_transfer[j];
    Edot_v_transfer[j] = Qdot_v_transfer[j];

end for;

/*** thermodynamic equilibrium at phase boundary ***/
for j in 1:n loop
   for i in 1:nS loop
    x_v_star_eq[j,mapping[i,1]]= K[j,i] *x_l_star[j,mapping[i,2]];
    //K[j,i] = thermoEquilibrium[j].K[i];
   end for;
end for;

if homotopyMethod.bool_K and homotopyMethod.useHomotopy then
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

  annotation (Documentation(info="<html>
<p><br/>There is the possiblity to assume thermodynamic equilibrium on each stage. This class sets the parameter EQ = true. If so, there is no resistance for mass transfer and the compositions in the bulk phase and at the phase boundary are the same. However for numerical reasons nevertheless a molar flow rate N_dot of each component is calculated using the equation</p>
<p>N_dot = A  k &middot; (x_star-x) for each component i.</p>
<p>The mass transfer coefficient k shall be chosen high enough to ensure that x_star &rarr; x for every component. However if K is too high the simulation becomes very slow or - even worse - the nonlinear solver fails to solve the problem. Also the optimal value for K may change during the simulation. Therefore only a start value for K is given by the user and then K is adapted using a PI controller. Up to now there is only one PI controller for the vapour phase and one for the liquid phase. The PI controller aims to get the maximum difference of the (x_star-x) of all substances on the first stage to zero.</p>
<p><br/><u><b><font style=\"color: #008000; \">Murhpree tray efficiency</font></b></u></p>
<p>This model also includes the Murphree tray efficiency. The Murphree efficiency &eta;Murphree is only used for tray columns (for packed columns the efficiency is set to one) and used together with an equilibrium model (parameter EQ=true): in this case it is said that the thermodynamic equilibrium is attained on every plate (plate column). However, if this is not the case but detailed information for a mass transfer model is missing, than the deviation from equilibrium can be expressed using the Murphree efficiency. The Murphree efficiency is defined as the ratio of the real concentration change on the tray to the maximal concentration change (which would be if equilibrium was obtained). The tray efficiency takes a concentration gradient in the liquid on the tray into account: The composition at the tray outlet is different than the composition on the tray (which is not modelled here: in the model for one stage the composition in the stage equals the composition at the stage outlet). Therefore it is possible for the tray efficiency to become greater than 1 (see for example Perry: Perry&apos;s chemical engineers&apos; handbook, 8th ed.). If no equilibrium model is used, &eta;Murphree is set to one. In a binary mixture &eta;Murphree is equal for both components; in a multicomponent mixture the efficiency is different for every component and can vary between -inf and inf.</p>
</html>"));
end EquilibriumEffectiveDiff;
