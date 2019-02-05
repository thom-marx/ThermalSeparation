within ThermalSeparation.BalanceEquations.Base.Equilibrium;
model BaseTwoPhaseSteadyState
  "Equilibrium: both phases balanced together, steady state"

extends
    ThermalSeparation.BalanceEquations.Base.Equilibrium.BaseBalanceEquationsEq;

  input SI.MolarInternalEnergy u_v[n](each stateSelect=StateSelect.default);
  input SI.MolarInternalEnergy u_l[n](each stateSelect=StateSelect.default);

  Modelica.SIunits.AmountOfSubstance n_mol_L[n](each stateSelect=StateSelect.default);
  Modelica.SIunits.AmountOfSubstance n_mol_V[n](each stateSelect=StateSelect.default);
  Modelica.SIunits.AmountOfSubstance n_mol[n,nS];

protected
  Real eps_liq_state[n](each stateSelect=StateSelect.default)=eps_liq;

equation
stat=true;

for j in 1:n loop
    n_mol_V[j] = A*H/n*eps*eps_vap[j]* sum(c_v[j,:]);
    n_mol_L[j] = A*H/n*eps*eps_liq[j]* sum(c_l[j,:]);
    for i in 1:nS loop
    n_mol[j,i]=A*H/n*eps*eps_vap[j]* c_v[j,mapping[i,1]]+A*H/n*eps*eps_liq[j]*c_l[j,mapping[i,1]];
    end for;
end for;

/*** MOLE BALANCES ***/
if n==1 then

/* EQ equations */

   // component balance for vapour and liquid together
  fill(0,nS) =  Ndot_v_in*x_upStreamIn_act + Ndot_l_in * x_downStreamIn_act  -Ndot_v[n]*x_upStreamOut_act -Ndot_l[1]* x_downStreamOut_act + Vdot_v_feed[1]*c_v_feed[1,:]+ Vdot_l_feed[1]*c_l_feed[1,:] + Ndot_source_startUp[1] * x_v[1,:]+ Ndot_reac[1,:]- Vdot_le[1]*c_l[1,:];

  // total mole balance for liquid and vapour together
   //bool_eps[1]=false;
  // der(sum(n_mol[1,:])) =  Vdot_l_in*rho_l_in/MM_l_in+Vdot_v_in*rho_v_in/MM_v_in -  Vdot_l[1]*rho_l[1]/MM_l[1]-  Vdot_v[1]*rho_v[1]/MM_v[1] + sum(Ndot_reac[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] + Vdot_v_feed[1]*rho_v_feed[1]/MM_v_feed[1]+ Ndot_source_startUp[1]- Vdot_le[1]*rho_l[1]/MM_l[1];
  sum(x_v[1,:])=1;
  sum(x_l[1,:])=1;
  Ndot_v_transfer[1,:]=fill(0,nS);

  else

    /** Begin lowest stage (n=1) **/
   // component balance for vapour and liquid together
   fill(0,nS) =  Ndot_v_in*x_upStreamIn_act + Ndot_l_in * x_downStreamIn_act  -Ndot_v[n]*x_upStreamOut_act -Ndot_l[1]* x_downStreamOut_act + Vdot_v_feed[1]*c_v_feed[1,:]+ Vdot_l_feed[1]*c_l_feed[1,:] + Ndot_source_startUp[1] * x_v[1,:]+ Ndot_reac[1,:]- Vdot_le[1]*c_l[1,:];

   // total mole balance for liquid and vapour together
   //bool_eps[1]=false;    /** End lowest stage (n=1) **/
   // der(sum(n_mol[1,:])) =  Vdot_l_in*rho_l_in/MM_l_in+Vdot_v_in*rho_v_in/MM_v_in -  Vdot_l[1]*rho_l[1]/MM_l[1]-  Vdot_v[1]*rho_v[1]/MM_v[1] + sum(Ndot_reac[1,:]) + Vdot_l_feed[1]*rho_l_feed[1]/MM_l_feed[1] + Vdot_v_feed[1]*rho_v_feed[1]/MM_v_feed[1]+ Ndot_source_startUp[1]- Vdot_le[1]*rho_l[1]/MM_l[1];
   sum(x_v[1,:])=1;
   sum(x_l[1,:])=1;
   Ndot_v_transfer[1,:]=fill(0,nS);

    /** Begin stages 2 to n-1 **/
    for j in 2:n-1 loop
      for i in 1:nSV loop
     0 = Vdot_v[j-1]*c_v[j-1,i]+Vdot_l[j+1]*c_l[j+1,i] - Vdot_v[j] * c_v[j,i]- Vdot_l[j] * c_l[j,i] + Ndot_reac[j,i]+ Vdot_l_feed[j]*c_l_feed[j,i] + Vdot_v_feed[j]*c_v_feed[j,i] + Ndot_source_startUp[j]*x_v[j,i]+ Vdot_le[j-1]*c_l[j-1,i] - Vdot_le[j]*c_l[j,i];
    end for;

   // total mole balance for liquid and vapour together
     //bool_eps[j]=false;
   // der(sum(n_mol[j,:])) = Vdot_l[j+1]*rho_l[j+1]/MM_l[j+1]+Vdot_v[j-1]*rho_v[j-1]/MM_v[j-1] -  Vdot_l[j]*rho_l[j]/MM_l[j]-  Vdot_v[j]*rho_v[j]/MM_v[j]+ sum(Ndot_reac[j,:]) + Vdot_l_feed[j]*rho_l_feed[j]/MM_l_feed[j]+ Vdot_v_feed[j]*rho_v_feed[j]/MM_v_feed[j] + Ndot_source_startUp[j] + Vdot_le[j-1]*rho_l[j-1]/MM_l[j-1] - Vdot_le[j]*rho_l[j]/MM_l[j];
   sum(x_v[j,:])=1;
   sum(x_l[j,:])=1;
   Ndot_v_transfer[j,:]=fill(0,nS);
    end for;
    /** End stages 2 to n-1 **/

    /** Begin highest stage (n=n) **/

   // component balance for vapour and liquid together
   fill(0,nS) = Vdot_v[n-1]*c_v[n-1,:]+Ndot_l_in * x_downStreamIn_act  -Ndot_v[n]*x_upStreamOut_act- Vdot_l[n] * c_l[n,:] + Vdot_v_feed[n]*c_v_feed[n,:]+ Ndot_reac[n,:] + Vdot_l_feed[n]*c_l_feed[n,:] + Vdot_le[n-1]*c_l[n-1,:]  + Ndot_source_startUp[n]*x_v[n,:];

  // total mole balance for liquid and vapour together
  // bool_eps[n]=false;
   // der(sum(n_mol[n,:])) = Ndot_l_in +Vdot_v[n-1]*rho_v[n-1]/MM_v[n-1] -  Vdot_l[n]*rho_l[n]/MM_l[n] -Ndot_v[n] + sum(Ndot_reac[n,:])+ Vdot_v_feed[n]*rho_v_feed[n]/MM_v_feed[n] + Ndot_source_startUp[n] + Vdot_l_feed[n]*rho_l_feed[n]/MM_l_feed[n] + Vdot_le[n-1]*rho_l[n-1]/MM_l[n-1];
   sum(x_v[n,:])=1;
   sum(x_l[n,:])=1;
   Ndot_v_transfer[n,:]=fill(0,nS);
  /** End highest stage (n=n) **/
  end if;

/*** ENERGY BALANCE ***/
  if n==1 then
    0 = Ndot_l_in*h_downStreamIn_act+Ndot_v_in*h_upStreamIn_act   -Ndot_l[1]*h_downStreamOut_act-Ndot_v[n]*h_upStreamOut_act  - Qdot_wall[1] + Vdot_l_feed[1]*sum(c_l_feed[1,:])*h_l_feed[1] - Vdot_le[1]*sum(c_l[1,:])*h_l[1] +Qdot_reac[1]+ Vdot_v_feed[1]*sum(c_v_feed[1,:])*h_v_feed[1] + Ndot_source_startUp[1]*h_v[1];
    Edot_l_transfer[1]=0;
  else
   0 = Vdot_l[2]*sum(c_l[2,:])*h_l[2] +Ndot_v_in*h_upStreamIn_act -Ndot_l[1]*h_downStreamOut_act  - Vdot_v[1]*sum(c_v[1,:])*h_v[1] - Qdot_wall[1] + Vdot_l_feed[1]*sum(c_l_feed[1,:])*h_l_feed[1] - Vdot_le[1]*sum(c_l[1,:])*h_l[1] + Qdot_reac[1]+ Vdot_v_feed[1]*sum(c_v_feed[1,:])*h_v_feed[1] + Ndot_source_startUp[1]*h_v[1];
   for j in 2:n-1 loop
     0 = Vdot_l[j+1]*sum(c_l[j+1,:])*h_l[j+1] +Vdot_v[j-1] * sum(c_v[j-1,:])*propsVap[j-1].h - Vdot_l[j]*sum(c_l[j,:])*h_l[j] - Vdot_v[j]*sum(c_v[j,:])*h_v[j] -Qdot_wall[j]+ Vdot_l_feed[j]*sum(c_l_feed[j,:])*h_l_feed[j] +Vdot_le[j-1]*sum(c_l[j-1,:])*h_l[j-1] - Vdot_le[j]*sum(c_l[j,:])*h_l[j] + Qdot_reac[j]+ Vdot_v_feed[j]*sum(c_v_feed[j,:])*h_v_feed[j]  + Ndot_source_startUp[j]*h_v[j];
   end for;
   0= Ndot_l_in*h_downStreamIn_act +Vdot_v[n-1] * sum(c_v[n-1,:])*propsVap[n-1].h - Vdot_l[n]*sum(c_l[n,:])*h_l[n] -Ndot_v[n]*h_upStreamOut_act - Qdot_wall[n] + Vdot_l_feed[n]*sum(c_l_feed[n,:])*h_l_feed[n] + Vdot_le[n-1]*sum(c_l[n-1,:])*h_l[n-1]+ Qdot_reac[n] + Vdot_v_feed[n]*sum(c_v_feed[n,:])*h_v_feed[n] + Ndot_source_startUp[n]*h_v[n];
   for i in 1:n loop
     Edot_l_transfer[i]=0;
   end for;

  end if;

end BaseTwoPhaseSteadyState;
