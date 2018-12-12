within ThermalSeparation.Components.ColumnsNoIndex.BaseClasses.Initialization;
model Init_p_Ndot0_Edot0_xtotal
  "initialization in equilibrium (Ndot_transfer = 0, Edot_transfer=0) using x_total and p"
extends BaseInit;
initial equation

  if considerStartUp then
        for j in 1:n loop
            for i in 1:nSL loop
           c_l[j,i]=  x_l_start[j,i] /v[j];
            end for;
            for i in 1:nSV loop
             x_v[j,i] = x_v_start[j,i];
            end for;
          end for;
          T_v = T_v_start;
          T_l = T_l_start;
        // p_hyd[1:n] = p_v_start;
         //oder
         p_v[1:n] = p_v_start;

else

    for j in 1:n loop
          for i in 1:nSL loop
             if inertLiquid[i] then
                x_l[j,i] = x_l_start[j,i];
             end if;
          end for;
       sum(x_l[j,:])=1;
       sum(x_v[j,:])=1;
    end for;
    p_v[1:n] = p_v_start;
    for i in 1:nSV loop
      if not inertVapour[i] then
        Ndot_v_transfer[:,i]=zeros(n);
      else
        x_v[:,i] = x_v_start[:,i];
      end if;
    end for;

         for j in 1:n loop
           for i in 1:nS-1 loop
             (n_mol_L[j]*x_l_star[j,mapping[i,2]]+n_mol_V[j]*x_v_star[j,mapping[i,1]])/(n_mol_L[j]+n_mol_V[j]) = x_total_start[mapping[i,1]];
           end for;
           end for;
Edot_l_transfer=zeros(n);

end if;
end Init_p_Ndot0_Edot0_xtotal;
