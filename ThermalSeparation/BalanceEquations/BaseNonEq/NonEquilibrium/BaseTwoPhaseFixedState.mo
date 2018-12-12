within ThermalSeparation.BalanceEquations.BaseNonEq.NonEquilibrium;
model BaseTwoPhaseFixedState
  "phases balanced seperately, states fixed, no index reduction performed"
extends ThermalSeparation.BalanceEquations.BaseNonEq.BaseBalanceEq;

/*** changes to reduce index: fixed state selection ***/
  output SI.MolarInternalEnergy u_v[n](stateSelect=StateSelect.always);
  output SI.MolarInternalEnergy u_l[n](stateSelect=StateSelect.always);
  Modelica.SIunits.AmountOfSubstance n_mol_L[n](stateSelect=StateSelect.always);
  Modelica.SIunits.AmountOfSubstance n_mol_V[n](stateSelect=StateSelect.always);
  Modelica.SIunits.AmountOfSubstance n_mol_L_i[n,nSL](stateSelect=StateSelect.always);
  Modelica.SIunits.AmountOfSubstance n_mol_V_i[n,nSV](stateSelect=StateSelect.always);
  Modelica.SIunits.AmountOfSubstance n_mol[n,nS];

equation
/*** MOLE BALANCES ***/
if n==1 then

   // component balance for vapour
  der(n_mol_V_i[1,:]) =  Ndot_v_in*x_upStreamIn_act  -Ndot_v[n]*x_upStreamOut_act + Ndot_v_transfer[1,:] + Vdot_v_feed[1]*c_v_feed[1,:] + Ndot_source_startUp[1] * x_v[1,:];
    // component balance for liquid
 der(n_mol_L_i[1,:]) = Ndot_l_in * x_downStreamIn_act  -Ndot_l[1]* x_downStreamOut_act + Ndot_l_transfer[1,:] + Ndot_reac[1,:] + Vdot_l_feed[1]*c_l_feed[1,:] - Vdot_le[1]*c_l[1,:];

  // total mole balance for liquid and vapour
   bool_eps[1]=false;
   der(n_mol_L[1]) =  Vdot_l_in*rho_l_in/MM_l_in -  Vdot_l[1]*rho_l[1]/MM_l[1] + sum(Ndot_l_transfer[1,:]) + sum(Ndot_reac[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];
   //    A*H/n*eps* der(eps_liq[1]*rho_l[1]/MM_l[1]) =  Vdot_l_in*rho_l_in/MM_l_in -  Vdot_l[1]*rho_l[1]/MM_l[1] + sum(Ndot_l_transfer[1,:]) + sum(reaction.Ndot[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];

      der(n_mol_V[1]) = Vdot_v_in*rho_v_in/MM_v_in -  Vdot_v[1]*rho_v[1]/MM_v[1] + sum(Ndot_v_transfer[1,:]) + Vdot_v_feed[1]*rho_v_feed[1]/MM_v_feed[1] + Ndot_source_startUp[1];
else

  /** Begin lowest stage (n=1) **/

   // component balance for vapour
  der(n_mol_V_i[1,:]) = Ndot_v_in*x_upStreamIn_act - Vdot_v[1] * c_v[1,:] + Ndot_v_transfer[1,:] + Vdot_v_feed[1]*c_v_feed[1,:] + Ndot_source_startUp[1]*x_v[1,:];
    // component balance for liquid
   der(n_mol_L_i[1,:]) = Vdot_l[2]*c_l[2,:]    -Ndot_l[1] * x_downStreamOut_act + Ndot_l_transfer[1,:] + Ndot_reac[1,:] + Vdot_l_feed[1]*c_l_feed[1,:] - Vdot_le[1]*c_l[1,:];

  // total mole balance for liquid and vapour

   bool_eps[1]=false;
   der(n_mol_L[1]) =  Vdot_l[2]*rho_l[2]/MM_l[2] -Ndot_l[1] + sum(Ndot_l_transfer[1,:]) + sum(Ndot_reac[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];
   //    A*H/n*eps* der(eps_liq[1]*rho_l[1]/MM_l[1]) =  Vdot_l[2]*rho_l[2]/MM_l[2] -  Vdot_l[1]*rho_l[1]/MM_l[1] + sum(Ndot_l_transfer[1,:]) + sum(reaction.Ndot[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] - Vdot_le[1]*rho_l[1]/MM_l[1];

   der(n_mol_V[1]) = Ndot_v_in -  Vdot_v[1]*rho_v[1]/MM_v[1] + sum(Ndot_v_transfer[1,:]) + Vdot_v_feed[1]*rho_v_feed[1]/MM_v_feed[1] + Ndot_source_startUp[1];
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

   der(n_mol_V_i[n,:]) = Vdot_v[n-1]*c_v[n-1,:]  -Ndot_v[n]*x_upStreamOut_act + Ndot_v_transfer[n,:] + Vdot_v_feed[n]*c_v_feed[n,:]  + Ndot_source_startUp[n]*x_v[n,:];

    // component balance for liquid
   der(n_mol_L_i[n,:]) = Ndot_l_in * x_downStreamIn_act - Vdot_l[n] * c_l[n,:] + Ndot_l_transfer[n,:]+ Ndot_reac[n,:] + Vdot_l_feed[n]*c_l_feed[n,:] + Vdot_le[n-1]*c_l[n-1,:];

  // total mole balance for liquid and vapour

     bool_eps[n]=false;
   der(n_mol_L[n]) = Ndot_l_in -  Vdot_l[n]*rho_l[n]/MM_l[n] + sum(Ndot_l_transfer[n,:])+ sum(Ndot_reac[n,:]) + Vdot_l_feed[n]*rho_l_feed[n]/MM_l_feed[n] + Vdot_le[n-1]*rho_l[n-1]/MM_l[n-1];
   //   A*H/n*eps* der(eps_liq[n]*rho_l[n]/MM_l[n]) = Vdot_l_in*rho_l_in/MM_l_in -  Vdot_l[n]*rho_l[n]/MM_l[n] + sum(Ndot_l_transfer[n,:])+ sum(reaction.Ndot[n,:]) + Vdot_l_feed[n]*rho_l_feed[n]/MM_l_feed[n] + Vdot_le[n-1]*rho_l[n-1]/MM_l[n-1];

   der(n_mol_V[n]) = Vdot_v[n-1]*rho_v[n-1]/MM_v[n-1] -Ndot_v[n] + sum(Ndot_v_transfer[n,:])+ Vdot_v_feed[n]*rho_v_feed[n]/MM_v_feed[n] + Ndot_source_startUp[n];
/** End highest stage (n=n) **/
end if;

/*** ENERGY BALANCE ***/
if n==1 then
    der(n_mol_L[1]*u_l[1]) =sum(-Ndot_v_transfer[1,:].* delta_hv[1,:])+ Ndot_l_in*h_downStreamIn_act   -Ndot_l[1]*h_downStreamOut_act  - Qdot_wall[1] + Edot_l_transfer[1]  + Vdot_l_feed[1]*sum(c_l_feed[1,:])*h_l_feed[1] - Vdot_le[1]*sum(c_l[1,:])*h_l[1] +Qdot_reac[1];
    der(n_mol_V[1]*u_v[1])  = Ndot_v_in*h_upStreamIn_act   -Ndot_v[n]*h_upStreamOut_act +Edot_v_transfer[1]  + Vdot_v_feed[1]*sum(c_v_feed[1,:])*h_v_feed[1] + Ndot_source_startUp[1]*h_v[1];
else
  der(n_mol_L[1]*u_l[1]) = sum(-Ndot_v_transfer[1,:].* delta_hv[1,:])+Vdot_l[2]*sum(c_l[2,:])*h_l[2] -Ndot_l[1]*h_downStreamOut_act - Qdot_wall[1] + Edot_l_transfer[1]  + Vdot_l_feed[1]*sum(c_l_feed[1,:])*h_l_feed[1] - Vdot_le[1]*sum(c_l[1,:])*h_l[1] + Qdot_reac[1];
    for j in 2:n-1 loop
  der(n_mol_L[j]*u_l[j]) = sum(-Ndot_v_transfer[j,:].* delta_hv[j,:])+  Vdot_l[j+1]*sum(c_l[j+1,:])*h_l[j+1] - Vdot_l[j]*sum(c_l[j,:])*h_l[j] -Qdot_wall[j]+Edot_l_transfer[j] + Vdot_l_feed[j]*sum(c_l_feed[j,:])*h_l_feed[j] +Vdot_le[j-1]*sum(c_l[j-1,:])*h_l[j-1] - Vdot_le[j]*sum(c_l[j,:])*h_l[j] + Qdot_reac[j];
   end for;
  der(n_mol_L[n]*u_l[n])=sum(-Ndot_v_transfer[n,:].* delta_hv[n,:])+ Ndot_l_in*h_downStreamIn_act - Vdot_l[n]*sum(c_l[n,:])*h_l[n] - Qdot_wall[n] + Edot_l_transfer[n] + Vdot_l_feed[n]*sum(c_l_feed[n,:])*h_l_feed[n] + Vdot_le[n-1]*sum(c_l[n-1,:])*h_l[n-1]+ Qdot_reac[n];

  der(n_mol_V[1]*u_v[1])  =  Ndot_v_in*h_upStreamIn_act - Vdot_v[1]*sum(c_v[1,:])*h_v[1] +Edot_v_transfer[1]  + Vdot_v_feed[1]*sum(c_v_feed[1,:])*h_v_feed[1] + Ndot_source_startUp[1]*h_v[1];
  for j in 2:n-1 loop
  der(n_mol_V[j]*u_v[j])  =Vdot_v[j-1] * sum(c_v[j-1,:])*propsVap[j-1].h  - Vdot_v[j]*sum(c_v[j,:])*h_v[j] +Edot_v_transfer[j] + Vdot_v_feed[j]*sum(c_v_feed[j,:])*h_v_feed[j]  + Ndot_source_startUp[j]*h_v[j];
  end for;
  der(n_mol_V[n]*u_v[n]) = Vdot_v[n-1] * sum(c_v[n-1,:])*propsVap[n-1].h  -Ndot_v[n]*h_upStreamOut_act +Edot_v_transfer[n] + Vdot_v_feed[n]*sum(c_v_feed[n,:])*h_v_feed[n] + Ndot_source_startUp[n]*h_v[n];
end if;

for j in 1:n loop
    n_mol_V[j] = A*H/n*eps*eps_vap[j]* sum(c_v[j,:]);
    n_mol_L[j] = A*H/n*eps*eps_liq[j]* sum(c_l[j,:]);
    n_mol_L_i[j,:]=n_mol_L[j]*x_l[j,:];
    n_mol_V_i[j,:]=n_mol_V[j]*x_v[j,:];
    for i in 1:nS loop
    n_mol[j,i]=A*H/n*eps*eps_vap[j]* c_v[j,mapping[i,1]]+A*H/n*eps*eps_liq[j]*c_l[j,mapping[i,1]];
    end for;
end for;

end BaseTwoPhaseFixedState;
