within ThermalSeparation.FilmModel;
package BaseClasses 
  model TrueEquilibriumStartUpCCSAbsorption "true equilibrium model for StartUp"
    extends BaseFilmModel(final EQ=true);

  Real K[n,nS] "equilibrium constant" annotation(Dialog(__Dymola_initialDialog=true));

    SI.Density rho_v[n]=propsVap.rho;
    SI.MolarMass MM_v[n] = propsVap.MM;

     replaceable model ThermoEquilibrium =
        PhaseEquilibrium.RealGasActivityCoeffLiquid constrainedby PhaseEquilibrium.BasePhaseEquilibrium
         annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

    ThermoEquilibrium thermoEquilibrium[n](
      each nS=nS,  each mapping =                                    mapping,
        redeclare replaceable package MediumVapour =  MediumVapour,
        redeclare replaceable package MediumLiquid =    MediumLiquid,
         p=p_v[1:n], T=T_star, x_v=x_v_star, x_l=x_l_star, p_sat=p_sat,  v_v=MM_v./rho_v,each x_vap_liq=fill(1/nS,nS));

    replaceable model StateSelection =
      ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.None
      constrainedby ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.BaseStateSelectionEq
                                                                                  annotation(choicesAllMatching);

    StateSelection stateSelection(
       redeclare replaceable package MediumVapour =  MediumVapour,
      redeclare replaceable package MediumLiquid =   MediumLiquid,
    n=n, propsLiq=propsLiq, propsVap=propsVap, c_l=c_l, p_v= p_v[1:n]);

    //        ThermalSeparation.Utilities.PI_Input PI_vap( initType=Modelica.Blocks.Types.Init.InitialOutput,
  //       T=T_vap,
  //       k=k_vap,
  //       y_start=K_start_vap, u= if startUp[1] then 0 else max(vector_vap2[1,:]));

    parameter Real k_vap = 1e1 "gain of PI controller" annotation(Dialog(tab="Heat and Mass Transfer", group = "Equilibrium"));
    parameter Real T_vap = 0.01 "time constant of PI controller" annotation(Dialog(tab="Heat and Mass Transfer", group = "Equilibrium"));
    parameter Real K_start_vap=2.5e-1 "start value for output of PI controller" annotation(Dialog(tab="Heat and Mass Transfer", group = "Equilibrium"));

   //Start-Up
    Real startUptime[n](start=fill(1e6,n));
    Real omega_smooth[n];
    parameter Real faktor_Ndot_v=300;

protected
     Real vector_vap2[ n,nSV];
     SI.MoleFraction x_v_star_eq[n,nSV];
  // Integer mapping2[nSV-nS];
  //Integer mapping2 [2]={2,4};
  //
  // Integer k;

  algorithm
  //   k:=1;
  //   for i in 1:nSV loop
  //     if inertVapour[i] then
  //   mapping2[k]:=i;
  //     else
  //       end if;
  //   k:=if inertVapour[i] then k + 1 else k;
  //   end for;
  equation
     for j in 1:n loop
    for i in 1:nSV loop
      vector_vap2[j,i] = (-x_v[j,i] + x_v_star[j,i]);
    end for;
    end for;

      /*** film is calculated steady-state ***/
    Ndot_v_interface=Ndot_v_transfer;
    Ndot_l_interface=Ndot_l_transfer;
    Edot_v_interface=Edot_v_transfer;
    Edot_l_interface=Edot_l_transfer;

        /*** energy balance at the phase boundary ***/
  for j in 1:n loop

  T_v[j] = T_star[j];
  T_l[j] = T_star[j];
  end for;

  if considerStartUp then
    for j in 1:n loop
       for i in 1:nSV loop
        if inertVapour[i] then
          x_v_star[j,i] = x_v[j,i];
        else
          if lowBoilingPoint[i] then
             x_v_star[j,i] = x_v[j,i];
          else
            if smooth_startUp then
                Ndot_v_transfer[j,i] = omega_smooth[j]*(-5.89)*faktor_Ndot_v*(x_v[j,i] - x_v_star[j,i]);
              else
                if  before_transition[j] then
                  Ndot_v_transfer[j,i] = 0;//if j==1 then -Vdot_v_in*c_v_in[i] else -Vdot_v[j-1]*c_v[j-1,i];//omega[j]*(  1.89*PI_vap.y*(x_v[j,i] - x_v_star[j,i])) + (1-omega[j])*(-Vdot_v_in*c_v_in[i]);// if j==1 then -Vdot_v_in*c_v_in[i] else -Vdot_v[j-1]*c_v[j-1,i];
                else
                  Ndot_v_transfer[j,i] = -1.89*faktor_Ndot_v*(x_v[j,i] - x_v_star[j,i]);
                end if;
            end if;
          end if;
        end if;
       end for;
    end for;
  else
     x_v_star = x_v;
  end if;

  /*** thermodynamic equilibrium at phase boundary ***/
  for j in 1:n loop
     for i in 1:nS loop
      x_v_star_eq[j,mapping[i,1]]= K[j,i] *x_l_star[j,mapping[i,2]];
     end for;
     for i in 1:nS loop
      K[j,i] = thermoEquilibrium[j].K[i];
     end for;
     for i in 1:nSV loop
       //x_v_star[j,i] = x_v[j,i];
     end for;
     for i in 1:nSL loop
       x_l_star[j,i] = x_l[j,i];
     end for;
  end for;

  for j in 1:n loop
    when startUp[j]==false then
      startUptime[j]= time;
    end when;
  end for;

  for j in 1:n loop
    if smooth_startUp then
      omega_smooth[j]=0.5 + 0.5*tanh(k*(time - (startUptime[j] + delay_startUp)));
    else
      omega_smooth[j]=0;
    end if;
  end for;

  //
  //  for j in 1:n loop
  //      if startUp[j]==false then
  //        assert(thermoEquilibrium[j].omega_init > 0.99, "Start-up condition is fulfilled while the initialization (homotopy) has not been completed yet");
  //      end if;
  //  end for;

    annotation (Documentation(info="<html>
<p>Thermodynamic equilibrium between the two bulk phases.</p>
</html>"));
  end TrueEquilibriumStartUpCCSAbsorption;

  model TrueEquilibriumStartUpCCSDesorption "true equilibrium model for Desorption StartUp"
    extends BaseFilmModel(final EQ=true);

  Real K[n,nS] "equilibrium constant" annotation(Dialog(__Dymola_initialDialog=true));

    SI.Density rho_v[n]=propsVap.rho;
    SI.MolarMass MM_v[n] = propsVap.MM;

     replaceable model ThermoEquilibrium =
        ThermalSeparation.PhaseEquilibrium.RealGasActivityCoeffLiquid constrainedby PhaseEquilibrium.BasePhaseEquilibrium
         annotation(Dialog(tab="Propagated from Column",group="These variables are propagated from the column model and do not have to be set by the user!",enable=false));

    ThermoEquilibrium thermoEquilibrium[n](
      each nS=nS,  each mapping =                                    mapping,
        redeclare replaceable package MediumVapour =  MediumVapour,
        redeclare replaceable package MediumLiquid =    MediumLiquid,
         p=p_v[1:n], T=T_star, x_v=x_v_star, x_l=x_l_star, p_sat=p_sat,  v_v=MM_v./rho_v, each x_vap_liq=fill(1/nS,nS));

    replaceable model StateSelection =
      ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.None
      constrainedby ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.BaseStateSelectionEq
                                                                                  annotation(choicesAllMatching);

    StateSelection stateSelection(
       redeclare replaceable package MediumVapour =  MediumVapour,
      redeclare replaceable package MediumLiquid =   MediumLiquid,
    n=n, propsLiq=propsLiq, propsVap=propsVap, c_l=c_l, p_v= p_v[1:n]);

    //        ThermalSeparation.Utilities.PI_Input PI_vap( initType=Modelica.Blocks.Types.Init.InitialOutput,
  //       T=T_vap,
  //       k=k_vap,
  //       y_start=K_start_vap, u= if startUp[1] then 0 else max(vector_vap2[1,:]));

    parameter Real k_vap = 1e1 "gain of PI controller" annotation(Dialog(tab="Heat and Mass Transfer", group = "Equilibrium"));
    parameter Real T_vap = 0.01 "time constant of PI controller" annotation(Dialog(tab="Heat and Mass Transfer", group = "Equilibrium"));
    parameter Real K_start_vap=2.5e-1 "start value for output of PI controller" annotation(Dialog(tab="Heat and Mass Transfer", group = "Equilibrium"));

    Real omega_smooth[n];
    Real startUptime[n](start=fill(1e6,n));
    parameter Real faktor_Ndot_v=100;

protected
     Real vector_vap2[ n,nSV];
     SI.MoleFraction x_v_star_eq[n,nSV];
  // Integer mapping2[nSV-nS];
  //Integer mapping2 [2]={2,4};
  //
  // Integer k;

  algorithm
  //   k:=1;
  //   for i in 1:nSV loop
  //     if inertVapour[i] then
  //   mapping2[k]:=i;
  //     else
  //       end if;
  //   k:=if inertVapour[i] then k + 1 else k;
  //   end for;
  equation
     for j in 1:n loop
    for i in 1:nSV loop
      vector_vap2[j,i] = (-x_v[j,i] + x_v_star[j,i]);
    end for;
    end for;

      /*** film is calculated steady-state ***/
    Ndot_v_interface=Ndot_v_transfer;
    Ndot_l_interface=Ndot_l_transfer;
    Edot_v_interface=Edot_v_transfer;
    Edot_l_interface=Edot_l_transfer;

   /*** energy balance at the phase boundary ***/
  for j in 1:n loop
    T_v[j] = T_star[j];
    T_l[j] = T_star[j];
  end for;

  for j in 1:n loop
    if smooth_startUp then
      omega_smooth[j]=0.5 + 0.5*tanh(k*(time - (startUptime[j] + delay_startUp)));
    else
      omega_smooth[j]=0;
    end if;
  end for;

  if considerStartUp and StartUp_CCS then
    for j in 1:n loop
       for i in 1:nSV loop
        if inertVapour[i] then
          x_v_star[j,i] = x_v[j,i];
        else
          if smooth_startUp then
            if lowBoilingPoint[i] then
              //Ndot_v_transfer[j,i] = omega_smooth[j]*(-5.89*300*(x_v[j,i] - x_v_star[j,i]));
              Ndot_v_transfer[j,i] = omega_smooth[j]*(-5.89*faktor_Ndot_v*(x_v[j,i] - x_v_star[j,i]));
            else
              //Ndot_v_transfer[j,i] = omega_smooth[j]*(-5.89*300*(x_v[j,i] - x_v_star[j,i])) + (1-omega_smooth[j])*(noEvent(if j==1 then -Vdot_v_in*c_v_in[i] else -Vdot_v[j-1]*c_v[j-1,i]));
              Ndot_v_transfer[j,i] = omega_smooth[j]*(-5.89*faktor_Ndot_v*(x_v[j,i] - x_v_star[j,i])) + (1-omega_smooth[j])*(noEvent(if j==1 then -Vdot_v_in*c_v_in[i] else -Vdot_v[j-1]*c_v[j-1,i]));
            end if;
          else
            if lowBoilingPoint[i] then
              if  startUp[j] then
                Ndot_v_transfer[j,i] = 0;//omega[j]*(  1.89*PI_vap.y*(x_v[j,i] - x_v_star[j,i])) + (1-omega[j])*(-Vdot_v_in*c_v_in[i]);// if j==1 then -Vdot_v_in*c_v_in[i] else -Vdot_v[j-1]*c_v[j-1,i];
              else
                Ndot_v_transfer[j,i] = -5.89*faktor_Ndot_v*(x_v[j,i] - x_v_star[j,i]);
              end if;
            else
              if  startUp[j] then
                Ndot_v_transfer[j,i] = if j==1 then -Vdot_v_in*c_v_in[i] else -Vdot_v[j-1]*c_v[j-1,i];//omega[j]*(  1.89*PI_vap.y*(x_v[j,i] - x_v_star[j,i])) + (1-omega[j])*(-Vdot_v_in*c_v_in[i]);// if j==1 then -Vdot_v_in*c_v_in[i] else -Vdot_v[j-1]*c_v[j-1,i];
              else
                Ndot_v_transfer[j,i] = -5.89*faktor_Ndot_v*(x_v[j,i] - x_v_star[j,i]);
              end if;
            end if;
          end if;
        end if;
       end for;
    end for;
  else
    x_v_star = x_v;
  end if;

  /*** thermodynamic equilibrium at phase boundary ***/
  for j in 1:n loop
     for i in 1:nS loop
      x_v_star_eq[j,mapping[i,1]]= K[j,i] *x_l_star[j,mapping[i,2]];
     end for;
     for i in 1:nS loop
      K[j,i] = thermoEquilibrium[j].K[i];
     end for;
     for i in 1:nSV loop
       //x_v_star[j,i] = x_v[j,i];
     end for;
     for i in 1:nSL loop
       x_l_star[j,i] = x_l[j,i];
     end for;
  end for;

  for j in 1:n loop
    when before_transition[j]==false then
      startUptime[j]= time;
    end when;
  end for;

  // for j in 1:n loop
  //     if startUp[j]==false then
  //       assert(thermoEquilibrium[j].omega_init > 0.99, "Start-up condition is fulfilled while the initialization (homotopy) has not been completed yet");
  //     end if;
  // end for;

    annotation (Documentation(info="<html>
<p>Thermodynamic equilibrium between the two bulk phases.</p>
</html>"));
  end TrueEquilibriumStartUpCCSDesorption;
end BaseClasses;
